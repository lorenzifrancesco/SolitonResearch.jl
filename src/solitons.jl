"""
Iterate the soliton finding routine over a set of equations.
Use precomputed values when possible
"""
function fill_solitons(;
  eqs = [GPE_1D, NPSE, NPSE_plus])

  sl = load_simulation_list()
  # sim_dict = [load_simulation("input/", GPE_1D)]
  # sim = sim_dict[1]
  # xspace!(sim.psi_0, sim)
  # kspace!(sim.psi_0, sim)
  # runsim(sim, info=true)
  # @info "survived"
  sl = filter(p -> (p.equation in eqs), sl)
  get_soliton.(sl, use_precomputed=true)
  return
end


"""
Compute the solitonic ground state for the given Simulation.
Compare results with ground state dictionary, and eventually use 
precomputed results. 
Also print the ground state in human readable CSV file"""
function get_soliton(
  sim::Sim;
  use_precomputed::Bool=true,
  info::Bool=false,
  save_path = "results/")

  sim.iswitch = -im
  eq = sim.equation
  gs_dict = load_soliton_dictionary()
  gamma_param = g2gamma(sim.g, sim.equation)
  time_requirement = @elapsed begin
    if haskey(gs_dict, hs(eq.name, gamma_param))
      if use_precomputed
        info && @info @sprintf("\t %8s:    |  x  ", sim.name)
      else
        delete!(gs_dict, hs(eq.name, gamma_param))
        sol = runsim(sim; info=info)
        @info "running"
        info && @info "total imaginary time $(sol.cnt * sim.dt)"
        push!(gs_dict, hs(eq.name, gamma_param) => sol.u)
      end
    else
      try
        sol = runsim(sim; info=info)
      catch err
        if isa(err, NpseCollapse)
          @warn "NPSE   Collapse"
        elseif isa(err, Gpe3DCollapse)
          @warn "3D-GPE Collapse"
        else
          throw(err)
        end
      end
      push!(gs_dict, hs(eq.name, gamma_param) => sol.u)
    end
  end
  info && @info @sprintf("Ground state time: %8.4fs", time_requirement)
  solution = gs_dict[hs(eq.name, gamma_param)]

  JLD2.save(join([save_path, "gs_dict.jld2"]), gs_dict)
  human_readable_soliton()
  nothing
end


"""
Load solitons from file and plot a single figure 
of superimposed ground states.
An passant, create a human readable version of the dictionary
"""
function plot_solitons(;
  save_path = "results/",
  media_path = "media/",
  save_plots::Bool=true,
  show_plots::Bool=false,
  info::Bool=false)

  soliton_dict = load_soliton_dictionary(info=info)
  human_readable_soliton()

  # pyplot(size=(359, 220))
  p = plot()
  @assert is_gamma_uniform(soliton_dict)
  gamma = ihs(first(soliton_dict)[1])[2]
  cnt = 1
  for (k, v) in soliton_dict
    # TODO improve the waste of time
    sl = load_simulation_list()
    sim = sl[cnt]
    solution::AbstractArray = v
    plot_final_density!(
      p,
      [solution],
      sim;
      label=sim.name,
      color=:black,
    )
    cnt+=1
  end
  @info "special case for gamma = 0.65"
  plot!(p, xlims=(-4, 4), ylims=(0.0, 0.6))
  plot!(
    p,
    xlabel=L"x",
    ylabel=L"|f|^2",
    legend=:topright,
    grid=false,
    smooth=true,
  )
  # display and save
  if show_plots
    display(p)
  end
  save_plots ? savefig(p, "media/" * string(gamma) * "_ground_states.pdf") :  nothing

  # zoomed version
  if show_plots
    display(p)
  end
  plot!(p, xlims=(-1, 1), ylims=(0.30, 0.52)) # 0.4, 0.45 for gamma 065
  save_plots ? savefig(p, "media/" * string(gamma) * "_ground_states_zoom.pdf") :  nothing
  nothing
end


"""
Check if all the solitons stored in the dictionary share the same gamma
"""
function is_gamma_uniform(soliton_dict)
  cnt = 1
  tmp = 0.0
  flag = true
  for (k, v) in soliton_dict
    @info ihs(k)
    if cnt == 1
      tmp = ihs(k)[2]
    end
    if tmp != ihs(k)[2]
      flag = false
    end
    cnt += 1
  end
  return flag
end


"""
Check if the dictionary exists, return it if it exists, 
  or save a new empty one and return it if doesn't
""" 
function load_soliton_dictionary(;
  save_path="results/",
  info::Bool=false
)
  if isfile(save_path * "gs_dict.jld2")
    gs_dict = JLD2.load(save_path * "gs_dict.jld2")
    info && @info "Found GS dictionary: " gs_dict
  else
    gs_dict = Dict{EquationType, Tuple{AbstractArray}}()
    JLD2.save(join([save_path, "gs_dict.jld2"]), gs_dict)
  end
  return gs_dict
end


"""
Simple method to utilize the ground states
"""
function get_ground_state(sim; info=false)
  sim.iswitch = -im
  res = Array(runsim(sim; info=info).u)
  @assert size(res) == sim.N
  return res
end

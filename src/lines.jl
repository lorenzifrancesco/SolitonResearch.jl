"""
  Set the simulations extremes in terms of velocity 
"""
function bar_interval(i)
  extremes = [0.0, 1.0]
  return extremes[i]
end

function vel_interval(i)
  extremes = [0.1, 1.0]
  return extremes[i]
end

"""
  Get the lines for every equation selected
"""
function fill_lines(
  gamma=0.65;
  use_precomputed_lines=false,
  eqs=[NPSE_plus])
  if Threads.nthreads() == 1
    @warn "running in single thread mode!"
  else
    @info "running in multi-thread mode: n_threads =" Threads.nthreads()
  end

  save_path = "results/"

  sl = load_simulation_list()
  sl = filter(p -> (p.equation in eqs), sl)
  prepare_for_collision!.(sl, gamma;
    use_precomputed_gs=true)

  # TODO to be inserted into function
  if isfile(save_path * "line_dict.jld2")
    @info "Loading Lines library..."
    line_dict = JLD2.load(save_path * "line_dict.jld2")
  else
    @info "No lines library found! Saving an empty one..."
    line_dict = Dict()
    JLD2.save(save_path * "line_dict.jld2", line_dict)
  end

  for sim in sl
    name = sim.equation.name
    @info "================Lining " name
    if haskey(line_dict, hs(name, gamma)) && use_precomputed_lines
      @info "Already found line for " name, gamma
    else
      # launch the line methods
      line = get_lines(sim, name;
        lines=1,
        sweep="bar",
        points=20)

      push!(line_dict, hs(name, gamma) => line)
      JLD2.save(save_path * "line_dict.jld2", line_dict)
    end
  end
end


function get_lines(
  sim::Sim{1,Array{Complex{Float64}}},
  name::String="noname";
  lines=2,
  sweep="vel",
  points=100,
  messages=true
)
  saveto = "../media/lines_$(name).pdf"
  max_vel = vel_interval(2)
  max_bar = bar_interval(2)
  @warn lines
  @warn points
  # asymmetric matrix: 
  @assert sweep in ["vel", "bar"]
  # TODO fix, very inefficient
  if lines > 1
    if sweep == "vel"
      vel_list = LinRange(vel_interval(1), max_vel, points)
      bar_list = LinRange(bar_interval(1), max_bar, lines)
      x_axis = vel_list
      y_axis = bar_list
    elseif sweep == "bar"
      vel_list = LinRange(vel_interval(1), max_vel, lines)
      bar_list = LinRange(bar_interval(1), max_bar, points)
      x_axis = bar_list
      y_axis = vel_list
    end
  elseif lines == 1
    if sweep == "vel"
      vel_list = LinRange(vel_interval(1), max_vel, points)
      bar_list = [bar_interval(1)]
      x_axis = vel_list
      y_axis = bar_list
    elseif sweep == "bar"
      vel_list = [vel_interval(1)]
      bar_list = LinRange(bar_interval(1), max_bar, points)
      x_axis = bar_list
      y_axis = vel_list
    end
  end
  tran = Array{Float64,2}(undef, (lines, points))
  refl = Array{Float64,2}(undef, (lines, points))

  @info "Filling sim grid..."
  sgrid = Array{Sim,2}(undef, (lines, points))
  archetype = sim

  # all sims have the same x
  mask_refl = map(xx -> xx > 0, archetype.X[1] |> real)
  mask_tran = map(xx -> xx < 0, archetype.X[1] |> real)
  avg_iteration_time = 0.0

  print("____________________________________________________________________\n")
  print("|tid|   x |   y |     dt|   T %|1-T-R %| collapse|   iter time|\n")
  print("____________________________________________________________________")
  full_time = @elapsed begin
    # TODO put threads in the inner loop
    Threads.@threads for ix in eachindex(x_axis)
      x = x_axis[ix]
      for (iy, y) in enumerate(y_axis)
        if sweep == "vel"
          loop_sim = imprint_vel_set_bar(archetype; vv=x, bb=y)
        elseif sweep == "bar"
          loop_sim = imprint_vel_set_bar(archetype; vv=y, bb=x)
        end
        collapse_occured = false
        sol = nothing
        maxim = -1.0
        this_iteration_time = 0.0
        try
          this_iteration_time = @elapsed (sol, maxim) = runsim(loop_sim; info=false)
          avg_iteration_time += this_iteration_time
        catch err
          if isa(err, NpseCollapse) || isa(err, Gpe3DCollapse)
            collapse_occured = true
          else
            throw(err)
          end
        end

        @assert loop_sim.manual == true
        if !collapse_occured
          final = sol.u[end]
          xspace!(final, loop_sim)
          tran[iy, ix] = ns(final, loop_sim, mask_tran)
          refl[iy, ix] = ns(final, loop_sim, mask_refl)
        else
          tran[iy, ix] = NaN
        end

        tile_mess = @sprintf("| %2i|  %3i|  %3i|  %.3f|",
          Threads.threadid(),
          ix,
          iy,
          loop_sim.dt
        )
        messages && print("\n" * tile_mess * @sprintf("   %3i|    %3i|      %s|%12.2f|",
                            collapse_occured ? 999 : Int(round(tran[iy, ix] * 100)),
                            collapse_occured ? 999 : Int(round((1 - tran[iy, ix] - refl[iy, ix]) * 100)),
                            collapse_occured ? "yes" : " no",
                            this_iteration_time
                          ))
      end
    end
    # TODO save lines
  end
  print("\n")
  @info "==================================================================="
  @info "Lines time    = " * @sprintf("%.3f", full_time)
  @info "% time in solver = " * @sprintf("%.3f, %.0f %% of lines time", avg_iteration_time, avg_iteration_time / full_time * 100)
  @info "Single tile time = " * @sprintf("%.3f", avg_iteration_time / lines * points)
  @info "==================================================================="
  print("\n")
  return tran
end


"""
in the 3D case we do not have sufficient GPU mem, so we go serially
"""
function get_lines(
  archetype::Sim{3,CuArray{Complex{Float64}}},
  name::String="noname";
  lines=2,
  sweep="vel",
  points=100,)
  saveto = "../media/lines_$(name).pdf"
  max_vel = vel_interval(2)
  max_bar = bar_interval(2)
  #
  # asymmetric matrix: 
  @assert sweep in ["vel", "bar"]
  if sweep == "vel"
    vel_list = LinRange(vel_interval(1), max_vel, points)
    bar_list = LinRange(bar_interval(1), max_bar, lines)
    # FIXME find a better way to do this 0.1->1.0
    x_axis = vel_list
    y_axis = bar_list
  elseif sweep == "bar"
    vel_list = LinRange(vel_interval(1), max_vel, lines)
    bar_list = LinRange(bar_interval(1), max_bar, points)
    x_axis = bar_list
    y_axis = vel_list
  end
  tran = Array{Float64,2}(undef, (lines, points))
  refl = Array{Float64,2}(undef, (lines, points))


  @info "Proceeding serially from the archetype..."
  # all sims have the same x
  mask_refl = map(xx -> xx > 0, archetype.X[1] |> real)
  mask_tran = map(xx -> xx < 0, archetype.X[1] |> real)

  @info "Running lining..."
  avg_iteration_time = 0.0
  iter = Iterators.product(enumerate(y_axis), enumerate(x_axis))
  full_time = @elapsed for ((iy, y), (ix, x)) in ProgressBar(iter)
    sim = deepcopy(archetype)
    collapse_occured = false
    if sweep == "vel"
      imprint_vel_set_bar!(sim; vv=x, bb=y)
    elseif sweep == "bar"
      imprint_vel_set_bar!(sim; vv=y, bb=x)
    end
    sol = nothing
    try
      avg_iteration_time += @elapsed sol = runsim(sim; info=false)
    catch err
      if isa(err, NpseCollapse) || isa(err, Gpe3DCollapse)
        collapse_occured = true
      else
        throw(err)
      end
    end
    # catch maxiters hit and set the transmission to zero
    if sim.manual == false
      if sol.retcode != ReturnCode.Success
        @info "Run complete, computing transmission..."
        @info "Detected solver failure"
        tran[iy, ix] = 0.0
        refl[iy, ix] = 1.0
        @info "T = " tran[iy, ix]
      else
        if !collapse_occured
          final = sol.u[end]
          @info "Run complete, computing transmission..."
          xspace!(final, sim)
          tran[iy, ix] = ns(final, sim, mask_tran)
          refl[iy, ix] = ns(final, sim, mask_refl)
        else
          @info "Run complete, detected collapse..."
          tran[iy, ix] = NaN
        end
        @info "T = " tran[iy, ix]
      end
    else
      if !collapse_occured
        final = sol.u[end]
        @info "Run complete, computing transmission..."
        xspace!(final, sim)
        tran[iy, ix] = ns(final, sim, mask_tran)
        refl[iy, ix] = ns(final, sim, mask_refl)
      else
        @info "Run complete, detected collapse..."
        tran[iy, ix] = NaN
        refl[iy, ix] = NaN
      end
      @info "T = " tran[iy, ix]
    end
    if !isapprox(tran[iy, ix] + refl[iy, ix], 1.0, atol=1e-5)
      @warn "T+R != 1.0"
    end
  end
  @info "Lining time            = " full_time
  @info "Total time in solver   = " avg_iteration_time
  @info "Average iteration time = " avg_iteration_time / lines^2

  JLD2.@save("tran_$(name).jld2", tran)
  JLD2.@save("refl_$(name).jld2", refl)
  return tran
end

"""
  Simple plotting of all the lines in the dictionary
"""
function view_all_lines(; sweep="bar")
  line_file = "results/line_dict.jld2"
  @assert isfile(line_file)
  ld = JLD2.load(line_file)
  for (k, v) in ld
    @info "found" ihs(k)
    @warn "check the size"
    if sweep == "vel"
      p = plot(xlabel="velocity", ylabel="T", title=ihs(k))
      x = LinRange(vel_interval(1), vel_interval(2), length(v[1, :]))
      y = LinRange(bar_interval(1), bar_interval(2), length(v[:, 1]))
    elseif sweep == "bar"
      p = plot(xlabel="barrier", ylabel="T", title=ihs(k))
      x = LinRange(bar_interval(1), bar_interval(2), length(v[1, :]))
      y = LinRange(vel_interval(1), vel_interval(2), length(v[:, 1]))
    end
    for iy = 1:size(v)[1]
      plot!(p, collect(x), v[iy, :], label=string(iy))
    end
    savefig(p, "media/lines_" * string(ihs(k)) * ".pdf")
    #  display(p)
  end
end


"""
  More sophisticated plotting
"""
function plot_all_lines(number_of_lines=4, sweep="bar")
  # pyplot(size=(350, 220))
  for i = 1:number_of_lines
    compare_all_lines(slice_choice=i, sweep=sweep)
  end
end

function compare_all_lines(; slice_choice=1, sweep="bar")
  # pyplot(size=(350, 220))
  line_file = "results/line_dict.jld2"
  @assert isfile(line_file)
  ld = JLD2.load(line_file)
  example = collect(values(ld))[1]
  if sweep == "vel"
    p = plot(xlabel=L"v", ylabel=L"T")
    x = LinRange(vel_interval(1), vel_interval(2), length(example[1, :]))
  elseif sweep == "bar"
    p = plot(xlabel=L"b", ylabel=L"T")
    x = LinRange(bar_interval(1), bar_interval(2), length(example[1, :]))
  end

  cnt = 1
  @info keys(ld)
  delete!(ld, hs("CQ", 0.65))
  delete!(ld, hs("G1", 0.65))
  for (k, v) in ld
    @info "found" ihs(k)
    # for iy in 1:size(v)[1]
    if true
      if slice_choice == 2
        if ihs(k)[1] == "G3"
          choice[32:end] = ones(50 - 31) * NaN
        elseif ihs(k)[1] == "Np"
          choice[39:end] = ones(50 - 38) * NaN
        end
      end
      if slice_choice == 1
        if ihs(k)[1] == "G3"
          choice[16:end] = ones(50 - 15) * NaN
        elseif ihs(k)[1] == "Np"
          choice[35:end] = ones(50 - 34) * NaN
        end
      end
    end
    plot!(
      p,
      collect(x),
      choice,
      linestyle=lineof(ihs(k)),
      color=colorof(ihs(k)),
      label=nameof(ihs(k)),
    )
    # end
    cnt += 1
  end
  plot!(p, grid=false, legend=:topright)
  savefig(p, "media/compare_lines_" * string(slice_choice) * ".pdf")
  return p
end


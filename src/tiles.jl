# JULIA_CUDA_SOFT_MEMORY_LIMIT = "95%"

"""
Iterate the tiles finding routine over a set of equations.
Use precomputed values when possible.

we can plot_finals values for debug, and return_maximum for 
collapse validation 
"""
function fill_tiles(;
  return_maximum=false,
  number_of_tiles=20,
  eqs=[GPE_1D],
  plot_finals=false,
  gamma=0.65
)
  if Threads.nthreads() == 1
    @warn "running in single thread mode!"
  end
  sl = load_simulation_list()
  sl = filter(p -> (p.equation in eqs), sl)
  get_tile.(sl,
    tiles=number_of_tiles,
    use_precomputed=false,
    return_maximum=return_maximum,
    plot_finals=plot_finals)
  plot_tiles()
end

"""
Get the tiles, check if they are precomputed 
(TODO get partially precomputed tile grid)
"""
function get_tile(
  sim::Sim{1,Array{Complex{Float64}}};
  use_precomputed=false,
  name::String="noname",
  return_maximum=false,
  tiles=100,
  plot_finals=false,
  messages=true,
  infos=false
)
  name = sim.equation.name
  gamma = g2gamma(sim.g, sim.equation)
  messages && @info "==============================================="
  messages && @info "\t\tTiling " * string(sim.equation.name)
  messages && @info @sprintf("\t\tUsing gamma: %.3f", gamma)
  messages && @info "==============================================="
  if plot_finals
    @warn "Plotting finals!"
  end
  saveto = "../media/tiles_$(name).pdf"
  max_vel = 1
  max_bar = 1
  #
  vel_list = LinRange(0.1, max_vel, tiles)
  bar_list = LinRange(0, max_bar, tiles)
  tran = Array{Float64,2}(undef, (tiles, tiles))
  refl = Array{Float64,2}(undef, (tiles, tiles))
  maxi = Array{Float64,2}(undef, (tiles, tiles))
  #initialize negative values
  tran = -0.1 * ones((tiles, tiles))
  refl = -0.1 * ones((tiles, tiles))
  maxi = -0.1 * ones((tiles, tiles))

  sim.iswitch = 1
  tile_dict = load_tile_dictionary()
  messages && @info "\tSetting ground state..."
  prepare_for_collision!(sim, gamma; use_precomputed_gs=true)
  if haskey(tile_dict, hs(name, gamma)) && use_precomputed
    messages && @info "Already found tile for " equation, gamma
  else
    @info "===============================================\n\t Running tiling subroutine"
    archetype = sim
    # all sims have the same x
    mask_refl = map(xx -> xx > 0, archetype.X[1] |> real)
    mask_tran = map(xx -> xx < 0, archetype.X[1] |> real)

    this_iteration_time = 0.0
    avg_iteration_time = 0.0
    counter = 0
    iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))

    print("____________________________________________________________________\n")
    print("|tid| num|   bx|   vx|     dt|   T %|1-T-R %| collapse|   iter time|\n")
    print("____________________________________________________________________")

    
    full_time = @elapsed begin
      Threads.@threads for vx in eachindex(vel_list)
        # @showprogress "Computing all the velocities..." for vx in eachindex(vel_list)
        vv = vel_list[vx]
        # messages && print("\t"*@sprintf("Free memory = %.3f GiB", Sys.free_memory() / 2^30))
        collapse_occured = false
        maxim = -1.0
        # @spawnat ipr+1 begin
        for (bx, bb) in enumerate(bar_list)
          loop_sim = imprint_vel_set_bar(archetype; vv=vv, bb=bb)
          if bx > 2 && isnan(tran[bx-1, vx]) && isnan(tran[bx-2, vx])
            # messages && @printf("\n Collapse shortcut!")
            collapse_occured = true
            this_iteration_time = 0.0
          end
          sol = nothing
          tile_mess = @sprintf("| %2i|%4i|  %3i|  %3i|  %.3f|",
            Threads.threadid(),
            bx + tiles * (vx - 1),
            vx,
            bx,
            loop_sim.dt
          )
          # messages && print("\n..."*tile_mess) 
          if !collapse_occured
            try
              this_iteration_time = @elapsed  (sol, maxim) = runsim(loop_sim; info=infos, return_maximum=return_maximum)
              avg_iteration_time += this_iteration_time
              # FIXME avoid NPSE+ memory filling problem
              GC.gc()

              if plot_finals
                pp = plot_final_density(
                  sol.u,
                  loop_sim;
                  show=false,
                  title=@sprintf("[vx=%3i, bx=%3i]/%3i", vx, bx, tiles)
                )
                savefig(pp, "media/checks/final_$(name)_$(vx)_$(bx)_$(tiles).pdf")
                qq = plot_axial_heatmap(
                  sol.u,
                  loop_sim.t,
                  loop_sim;
                  show=false,
                  title=@sprintf("[vx=%i, bx=%i]/%i", vx, bx, tiles)
                )
                savefig(qq, "media/checks/heatmap_$(name)_$(vx)_$(bx)_$(tiles).pdf")
              end
            catch err
              if isa(err, NpseCollapse) || isa(err, Gpe3DCollapse)
                collapse_occured = true
              else
                throw(err)
              end
            end
          end
          @assert loop_sim.manual == true
          if !collapse_occured
            final = sol.u[end]
            # @info "Run complete, computing transmission..."
            xspace!(final, loop_sim)
            tran[bx, vx] = ns(final, loop_sim, mask_tran)
            refl[bx, vx] = ns(final, loop_sim, mask_refl)
            maxi[bx, vx] = maxim
          else
            tran[bx, vx] = NaN
            maxi[bx, vx] = NaN
          end
          messages && print("\n" * tile_mess * @sprintf("   %3i|    %3i|      %s|%12.2f|",
                              collapse_occured ? 999 : Int(round(tran[bx, vx] * 100)),
                              collapse_occured ? 999 : Int(round((1 - tran[bx, vx] - refl[bx, vx]) * 100)),
                              collapse_occured ? "yes" : " no",
                              this_iteration_time
                            ))
          counter += 1

        end # barrier loop
      # end # spawnat
      # ipr += 1
      # ipr = ipr % workers()

      end # velocities loop
      CSV.write("results/runtime_tran3.csv", Tables.table(tran))
      # csv2color("runtime_tran")
      if return_maximum
        CSV.write("results/runtime_maxi3.csv", Tables.table(maxi))
        # csv2color("runtime_maxi")
      end
    end
    messages && @info "Saving tiles"
    push!(tile_dict, hs(name, gamma) => tran)
    JLD2.save("results/tile_dict.jld2", tile_dict)
    print("\n")
    @info "==================================================================="
    @info "Pavement time    = " * @sprintf("%.3f", full_time)
    @info "% time in solver = " * @sprintf("%.3f, %.0f %% of pavement time", avg_iteration_time, avg_iteration_time / full_time * 100)
    @info "Single tile time = " * @sprintf("%.3f", avg_iteration_time / tiles^2)
    @info "==================================================================="
    print("\n")
  end
  return tran
end


"""
in the 3D case we do not have sufficient GPU mem, so we go serially
"""
function get_tile(
  archetype::Sim{3,CuArray{Complex{Float64}}},
  name::String="noname";
  tiles=100,
  plot_finals=false,
)
  @assert false
  saveto = "../media/tiles_$(name).pdf"
  max_vel = 1
  max_bar = 1
  #
  vel_list = LinRange(0.1, max_vel, tiles)
  bar_list = LinRange(0, max_bar, tiles)
  tran = Array{Float64,2}(undef, (tiles, tiles))
  refl = Array{Float64,2}(undef, (tiles, tiles))

  @info "Proceeding serially from the archetype..."
  # all sims have the same x
  mask_refl = map(xx -> xx > 0, archetype.X[1] |> real)
  mask_tran = map(xx -> xx < 0, archetype.X[1] |> real)

  @info "Running tiling..."
  avg_iteration_time = 0.0
  iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))
  full_time = @elapsed for ((vx, vv), (bx, bb)) in ProgressBar(iter)
    sim = deepcopy(archetype)
    collapse_occured = false
    imprint_vel_set_bar!(sim; vv=vv, bb=bb)
    @info "Computing tile" (vv, bb)
    sol = nothing
    try
      avg_iteration_time += @elapsed sol = runsim(sim; info=false)
      if plot_finals
        pp = plot_final_density(sol.u, sim; show=false)
        savefig(pp, "media/checks/final_$(name)_$(vv)_$(bb).pdf")
        qq = plot_axial_heatmap(sol.u, sim.t, sim; show=false)
        savefig(qq, "media/checks/heatmap_$(name)_$(vv)_$(bb).pdf")
      end

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
        tran[bx, vx] = 0.0
        refl[bx, vx] = 1.0
        @info "T = " tran[bx, vx]
      else
        if !collapse_occured
          final = sol.u[end]
          @info "Run complete, computing transmission..."
          xspace!(final, sim)
          tran[bx, vx] = ns(final, sim, mask_tran)
          refl[bx, vx] = ns(final, sim, mask_refl)
        else
          @info "Run complete, detected collapse..."
          tran[bx, vx] = NaN
        end
        @info "T = " tran[bx, vx]
      end
    else
      if !collapse_occured
        final = sol.u[end]
        @info "Run complete, computing transmission..."
        xspace!(final, sim)
        tran[bx, vx] = ns(final, sim, mask_tran)
        refl[bx, vx] = ns(final, sim, mask_refl)
      else
        @info "Run complete, detected collapse..."
        tran[bx, vx] = NaN
        refl[bx, vx] = NaN
      end
      @info "T = " tran[bx, vx]
    end
    if !isapprox(tran[bx, vx] + refl[bx, vx], 1.0, atol=1e-5)
      @warn "T+R != 1.0"
    end
  end
  @info "Tiling time            = " full_time
  @info "Total time in solver   = " avg_iteration_time
  @info "Average iteration time = " avg_iteration_time / tiles^2

  JLD2.@save("tran_$(name).jld2", tran)
  JLD2.@save("refl_$(name).jld2", refl)
  norm_bar = bar_list / max_bar
  norm_vel = vel_list / max_vel
  return tran
  return tran
end


"""
Check if the dictionary exists, return it if it exists, 
  or save a new empty one and return it if doesn't
"""
function load_tile_dictionary(;
  save_path="results/",
  info::Bool=false
)
  if isfile(save_path * "tile_dict.jld2")
    tile_dict = JLD2.load(save_path * "tile_dict.jld2")
    info && @info "Found tile dictionary: " tile_dict
  else
    tile_dict = Dict{EquationType,Tuple{AbstractArray}}()
    JLD2.save(join([save_path, "tile_dict.jld2"]), tile_dict)
  end
  return tile_dict
end
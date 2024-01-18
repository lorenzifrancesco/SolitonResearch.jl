using SolitonResearch, SolitonDynamics
using Plots, Printf, DataFrames, CSV
gr()
  
function precision(vx, bx, n_tiles=50) 
  N_samples = 1
  if N_samples > 1
    dt_list = exp(1) .^ LinRange(log(0.1), log(0.05), N_samples)
  else
    dt_list = [0.01]
  end
  vv = vx/n_tiles
  bb = bx/n_tiles
  inner_product = zeros(N_samples)
  exec_time = zeros(N_samples)

  sl = load_simulation_list()
  sl = filter(p -> (p.equation in [NPSE_plus]), sl)
  sim = sl[1]
  prepare_for_collision!(sim, 0.65, use_precomputed_gs=true)
  unew = zeros(length(sim.psi_0))
  uold = zeros(length(sim.psi_0))
  for (id, dd) in enumerate(dt_list)
    try
      exec_time[id] = @elapsed unew = get_final(sim, vv, bb, dd)
      if id == 1
        uold = unew
      end
      inner_product[id] = distance(unew, uold)
      uold = unew
    catch err 
      throw(err)
      @warn "collapse"
      inner_product[id] = -1
    end
    df = DataFrame([[dt_list]; [inner_product]; [exec_time]], :auto)
    display(df)
  end
  q = plot(dt_list, inner_product)
  savefig(q, "diagnostic/media/convergence.pdf")
  df = DataFrame([[dt_list]; [inner_product]; [exec_time]], :auto)
  display(df)
  CSV.write("diagnostic/dt_iterations.csv", df)
end

function get_final(sim, vx, bx, n_tiles, dt_set)
    ## load the simulations
    # sd = load_parameters()
    ## prepare in gs
    # prepare_for_collision!(sd, 0.65, use_precomputed_gs=true)
    ## select the NPSE+
    ## imprint velocity set barrier
  
    # vel = 0.1
    # if sim.equation == NPSE_plus
    #   x0 = sim.L[1]/8
    # else
    #   x0 = sim.L[1]/4
    # end
    # tf = 2*x0/vel
    # N = 50
    # dt_N = tf/N
    # @warn tf
    # @warn dt_N
    # @warn tf/dt_N
    imprint_vel_set_bar!(sim,
      dt_set=dt_set,
      vv=vx/n_tiles,
      bb=bx/n_tiles,
      save_each = true
    )
    # interesting case (0.3, 0.11111)
    ## print some infos and run the simulations
    # @info sim.dt
    # @info sim.tf
    # @info sim.L[1]
    # show_psi_0(sim)  
    sol = runsim(sim, info=true)

    psi2 = abs2.(xspace(sol.u[end], sim))
    p = plot(real(sim.X[1]), psi2)
    savefig(p, "diagnostic/media/tmp_vx_"* string(vx)*"_bx_"*string(bx)*"_final.pdf")
    plot_axial_heatmap(sol.u, sol.t, sim, path="diagnostic/media/")
    sleep(0.2)
  return sol.u[end]
end

function distance(u1, u2)
  @assert length(u1)==length(u2)
  return sum(abs2.(u1 .- u2))/sum(abs2.(u1))
end

precision()
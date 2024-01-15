using SolitonResearch, SolitonDynamics
using Plots, Printf
gr()

begin
  ## load the simulations
  sd = load_parameters()
  ## prepare in gs
  prepare_for_collision!(sd, 0.65, use_precomputed_gs=true)
  ## select the NPSE+
  sim = sd["Np"]
  ## imprint velocity set barrier
  imprint_vel_set_bar!(sim,
    dt_set=0.01,
    vv=0.4,
    bb=0.33333333
  )
  # interesting case (0.3, 0.11111)

  ## print some infos and run the simulations
  @info sim.dt
  @info sim.tf
  @info sim.L[1]
  # show_psi_0(sim)
end

@time begin
  sol = runsim(sim, info=true)
end

begin
  psi2 = abs2.(xspace(sol.u[end], sim))
  p = plot(real(sim.X[1]), psi2)
  savefig(p, "diagnostic/media/tmp_final.pdf")
  plot_axial_heatmap(sol.u, sol.t, sim, path="diagnostic/media/")
end
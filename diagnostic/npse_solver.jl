using SolitonResearch, SolitonDynamics
using Plots; pyplot()

function runna(sim, vx, bx)
  ## select the nasty simulation parameters
  total = 20
  imprint_vel_set_bar!(sim, vv=vx/total, bb=bx/total)
  # sim.time_steps=100
  return runsim(sim, info=true, return_maximum=true)
end


function lesgo(
  vx=20,
  bx=15
  )
  sim = load_simulation("input/", NPSE_plus)
  sim.Nt
  prepare_for_collision!(sim, 0.65, use_precomputed_gs=true)
  # show_psi_0(sim)
  @time (sol, maxi) = runna(sim, vx, bx)
  @info "======== ENDED ========="

  density_mia = abs2.(xspace(sol.u[end], sim))
  p = plot(real(sim.X[1]), density_mia)
  pht = plot_axial_heatmap(sol.u, sol.t, sim)
  show(p)
  q = plot(real(sim.X[1]), sol.sigma[end])
  show(q)
  qht = plot_axial_heatmap(sol.sigma, sol.t, sim; doifft=false)
  savefig(pht, "PSI__"*string(vx)*"_"*string(bx)*".pdf")
  savefig(qht, "SIGMA"*string(vx)*"_"*string(bx)*".pdf")
  return p, q
end
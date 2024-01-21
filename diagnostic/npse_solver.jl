using SolitonResearch, SolitonDynamics
using Plots

function runna(sim, vx, bx)
  ## select the nasty simulation parameters
  total = 20
  imprint_vel_set_bar!(sim, vv=vx/total, bb=bx/total)
  return runsim(sim, info=true)
end


function lesgo(
  vx=20,
  bx=4
  )
  sim = load_simulation("input/", NPSE_plus)
  sim.Nt
  prepare_for_collision!(sim, 0.65)
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
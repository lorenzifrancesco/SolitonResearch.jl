using SolitonResearch, SolitonDynamics
using Plots

function runna(vx, bx)
  ## select the nasty simulation parameters
  total = 20
  imprint_vel_set_bar!(sim, vv=vx/total, bb=bx/total)
  return runsim(sim, info=true)
end

sim = load_simulation("input/", NPSE_plus)
sim.Nt
prepare_for_collision!(sim, 0.65)

#####
@time (sol, maxi) = runna(20, 7)
density_mia = abs2.(xspace(sol.u[end], sim))
p = plot(real(sim.X[1]), density_mia)
plot_axial_heatmap(sol.u, sol.t, sim)
show(p)
q = plot(real(sim.X[1]), sol.sigma[end])
show(q)
plot_axial_heatmap(sol.sigma, sol.t, sim; doifft=false)
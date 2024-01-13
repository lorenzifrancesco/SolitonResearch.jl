using SolitonResearch, SolitonDynamics
using Plots, Printf
gr()

## load the simulations
sd = load_parameters()
## prepare in gs
prepare_for_collision!(sd, 0.65)
## select the NPSE+
sim = sd["Np"]
## imprint velocity set barrier
imprint_vel_set_bar!(sim,
  dt_set=0.01,
  vv=0.4,
  bb=0.333333333
)

## print some infos and run the simulations
@info sim.dt
@info sim.tf
@info sim.L[1]

@time begin
  sol = runsim(sim)
end

psi2 = abs2.(xspace(sol.u[end], sim))
p = plot(real(sim.X[1]), psi2)
savefig(p, "diagnostic/media/tmp.pdf")

plot_axial_heatmap(sol.u, sol.t, sim)
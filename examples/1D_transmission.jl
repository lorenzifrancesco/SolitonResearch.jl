import GR
using CondensateDynamics, CUDA, FFTW, Makie, GLMakie
using OrdinaryDiffEq
using LSODA
import CondensateDynamics.V

include("plot_axial_evolution.jl")

L = (70.0,)
N = (1024,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)


# ======= packing sim
@unpack_Sim sim
# ======= simulation custom parameters
equation = NPSE_plus
manual = true
time_steps = 500
Nt = 200

solver = SplitStep
g = -1.17
gamma_damp = 0.0
tiles = 8
barrier_width = 0.699 # as in SolitonBEC.jl
max_vel = 1.17 # CALCULATED VALUE 1.17 FOR CHOOSEN NONLINEARITY
max_bar = 1.68 # CALCULATED VALUE 1.68 FOR CHOOSEN NONLINEARITY
vx = 8
bx = 8

reltol = 1e-4
abstol = 1e-4

x0 = L[1] / 4

# other computations
iswitch = 1
g_param = abs(g) /2
sigma2 = init_sigma2(g)
x = X[1] |> real
k = K[1] |> real
dV= volume_element(L, N)

vel_list = LinRange(0, max_vel, tiles)
bar_list = LinRange(0, max_bar, tiles)
vv = vel_list[vx]
bb = bar_list[bx]
if vv == 0.0
    tf = 2.0
else
    tf = 2*x0/vv
end
t = LinRange(ti, tf, Nt)

dt = (tf-ti)/time_steps
@info "dt" dt
maxiters = 20000

@. psi_0 = sqrt(g_param/2) * 2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param)) * exp(-im*(x)*vv)

psi_0 = psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = BS3()

@. V0 = bb*exp(-(x/barrier_width)^2)

@pack_Sim! sim
@info sim
mask_tran = map(xx -> xx<0, x)

sol = runsim(sim; info=false)
if isnothing(sol)
    throw("NPSE collapse detected, cannot proceed further to plots...")
end

final = sol.u[end]
final = xspace(final, sim)
tran = ns(final, sim, mask_tran)
@info "T = " tran

time_axis = sol.t |> real
@info "Plotting..."
#plot_final_density(sol.u, psi_0, sim)
plot_axial_heatmap(sol.u, time_axis, sim)

# sigma heatmaps
plot_axial_heatmap(sigma2_new , time_of_sigma, sim; doifft = false)
plot_axial_heatmap(sigma2_old, time_of_sigma, sim; doifft = false)

# @info "Building animation..."
# animation_final_density(sol.u, sim)
# animation_final_density(sigma2_new, sim; doifft=false, info=true, file="sigma2_new.gif")
# animation_final_density(sigma2_old, sim; doifft=false, info=true, file="sigma2_old.gif")
# @info "Done!"

display(sol.destats)
display(sol.retcode)
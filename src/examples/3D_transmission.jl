using CondensateDynamics, OrdinaryDiffEq, LSODA
import CondensateDynamics.V
using CUDA
using LaTeXStrings, Plots
import GR
using CUDA.CUFFT
import Makie, GLMakie

include("plot_isosurfaces.jl")
update_parameters()
include("plot_axial_evolution.jl")

use_precomputed = false

# =================== simulation settings
L = (40.0,40.0,40.0)
N = (128,128,128)
sim = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)

# =================== physical parameters 
@unpack_Sim sim
g = -1.1 # this is g_1D (g_3D is computed in the propagation routine)
g_param = abs(g)/2
equation = GPE_3D
manual = true
iswitch = 1
x = Array(X[1])
y = Array(X[2])
z = Array(X[3])
dV= volume_element(L, N)
reltol = 1e-3
time_steps = 50
Nt = 30

# nfiles = true
maxiters = 20000
tiles = 8
barrier_width = 0.699 # as in SolitonBEC.jl
max_vel = 1.17 # CALCULATED VALUE 1.17 FOR CHOOSEN NONLINEARITY
max_bar = 1.68 # CALCULATED VALUE 1.68 FOR CHOOSEN NONLINEARITY
vx = 8  # 6
bx = 8 #  3
x0 = L[3]/4

vel_list = LinRange(0, max_vel, tiles)
bar_list = LinRange(0, max_bar, tiles)
vv = vel_list[vx]
bb = bar_list[bx]
tf = x0*2/vv
t = LinRange(ti,tf,Nt)
dt = (tf-ti)/time_steps

tmp = [exp(-(x^2+y^2)/2) * sqrt(g_param/2) * 2/(exp(g_param*(z-x0)) + exp(-(z-x0)*g_param)) * exp(-im*(z-x0)*vv) for x in x, y in y, z in z]
psi_0 = CuArray(tmp)

psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = BS3()
tmp = [1/2*(x^2+y^2) + bb*exp(-(z/barrier_width)^2) for x in x, y in y, z in z]
V0 = CuArray(tmp)
#V(x,y,z,t) = 1/2 * (x^2+y^2+z^2)
@pack_Sim! sim


# ===================== simulation
@info "computing GPE_3D" 
if isfile(join([save_path, "3d_tran.jld2"])) && use_precomputed
    @info "\t using precomputed solution 3d_tran.jld2" 
    JLD2.@load join([save_path, "3d_tran.jld2"]) gpe_1d
else
    sol = runsim(sim; info=false)
    u = sol.u
    JLD2.@save join([save_path, "3d_tran.jld2"]) u
end

# =================== plotting and collect 
#plot_final_density(u, psi_0, sim, 3; info=true)
plot_axial_heatmap(u, sol.t, sim, 3)

@info "Building animation..."
isosurface_animation(sol.u,length(sol.u), sim; framerate=5)
@info "Completed."
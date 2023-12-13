import GR
using CondensateDynamics
using OrdinaryDiffEq
using LSODA
import CondensateDynamics.V

include("plot_axial_evolution.jl")
save_path = "../results/"

gamma_param = 0.6
use_precomputed = false
maxiters_1d = 1000
N_axial_steps = 512

Nt = 30 # number of saves


####
tiles = 10
barrier_width = 0.699 # as in SolitonBEC.jl
max_vel = 1.17 # CALCULATED VALUE 1.17 FOR CHOOSEN NONLINEARITY
max_bar = 1.68 # CALCULATED VALUE 1.68 FOR CHOOSEN NONLINEARITY
vel_list = LinRange(0, max_vel, tiles)
bar_list = LinRange(0, max_bar, tiles)
vx = 8
bx = 8
####


# =========================================================
## 1D-GPE 

L = (40.0,)
N = (N_axial_steps,)
sim_gpe_1d = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim_gpe_1d
iswitch = 1
equation = GPE_1D
manual = true
solver = SplitStep
x0 = L[1] / 4
Nt = 30 # number of saves

# interaction parameter
g_param = gamma_param
maxiters = maxiters_1d
g = - 2 * g_param

abstol = 1e-6
x = X[1]
k = K[1]
dV= volume_element(L, N)
flags = FFTW.EXHAUSTIVE

vv = vel_list[vx]
bb = bar_list[bx]
if vv == 0.0
    tf = 2.0
else
    tf = 2*x0/vv
end
t = LinRange(ti, tf, Nt)
@. psi_0 = sqrt(g_param/2) * 2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param)) * exp(-im*(x-x0)*vv)
@. V0 = bb*exp(-(x/barrier_width)^2)

kspace!(psi_0, sim_gpe_1d)

@pack_Sim! sim_gpe_1d

# =========================================================
## NPSE (unable to copy)
L = (40.0,)
N = (N_axial_steps,)
sim_npse = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim_npse
iswitch = 1
equation = NPSE
manual = false
solver = SplitStep

# interaction parameter
g_param = gamma_param
maxiters = maxiters_1d
g = - 2 * g_param

abstol = 1e-6
x = X[1]
k = K[1]
dV= volume_element(L, N)
flags = FFTW.EXHAUSTIVE

vv = vel_list[vx]
bb = bar_list[bx]
if vv == 0.0
    tf = 2.0
else
    tf = 2*x0/vv
end
t = LinRange(ti, tf, Nt)
@. psi_0 = sqrt(g_param/2) * 2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param)) * exp(-im*(x-x0)*vv)
@. V0 = bb*exp(-(x/barrier_width)^2)
kspace!(psi_0, sim_gpe_1d)
if g_param > 2/3
    @warn "we should expect NPSE collapse"
end
sigma2 = init_sigma2(g)
@pack_Sim! sim_npse

# =========================================================
## NPSE (unable to copy)
L = (40.0,)
N = (N_axial_steps,)
sim_npse_plus = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim_npse_plus
iswitch = 1
equation = NPSE_plus
manual = false
solver = SplitStep
time_steps = 0

# interaction parameter
g_param = gamma_param
maxiters = maxiters_1d
g = - 2 * g_param

abstol = 1e-6
dt = 0.022
x = X[1]
k = K[1]
dV= volume_element(L, N)
flags = FFTW.EXHAUSTIVE

vv = vel_list[vx]
bb = bar_list[bx]
if vv == 0.0
    tf = 2.0
else
    tf = 2*x0/vv
end
t = LinRange(ti, tf, Nt)
@. psi_0 = sqrt(g_param/2) * 2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param)) * exp(-im*(x-x0)*vv)
@. V0 = bb*exp(-(x/barrier_width)^2)

kspace!(psi_0, sim_gpe_1d)
if g_param > 2/3
    @warn "we should expect NPSE collapse"
end
sigma2 = init_sigma2(g)
@pack_Sim! sim_npse_plus

# =========================================================
## 3D-GPE 

L = (40.0,40.0,40.0)
N = (N_axial_steps, 128, 128)
sim_gpe_3d = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim_gpe_3d
iswitch = 1
equation = GPE_3D
manual = false
solver = SplitStep

time_steps = 10000
g = - g_param * (4*pi)

reltol = 1e-3
maxiters = 500
alg = BS3()

x = Array(X[1])
y = Array(X[2])
z = Array(X[3])
dV= volume_element(L, N)

flags = FFTW.EXHAUSTIVE
vv = vel_list[vx]
bb = bar_list[bx]
if vv == 0.0
    tf = 2.0
else
    tf = 2*x0/vv
end
t = LinRange(ti, tf, Nt)

tmp = [exp(-(x^2+y^2)/2) * sqrt(g_param/2) * 2/(exp(g_param*(z-x0)) + exp(-(z-x0)*g_param)) * exp(-im*(z-x0)*vv) for x in x, y in y, z in z]
psi_0 = CuArray(tmp)
psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))

kspace!(psi_0, sim)

tmp = [1/2*(x^2+y^2) + bb*exp(-(z/barrier_width)^2) for x in x, y in y, z in z]
V0 = CuArray(tmp)

@pack_Sim! sim_gpe_3d

# =========================================================
time_axis = t
Plots.CURRENT_PLOT.nullableplot = nothing

@info "computing GPE_1D" 
if isfile(join([save_path, "gpe_1d_heat.jld2"])) && use_precomputed
    @info "\t using precomputed solution gpe_1d_heat.jld2" 
    JLD2.@load join([save_path, "gpe_1d_heat.jld2"]) gpe_1d
else
    sol = runsim(sim_gpe_1d; info=false)
    gpe_1d = sol.u
    JLD2.@save join([save_path, "gpe_1d_heat.jld2"]) gpe_1d
end
p = plot_axial_heatmap(gpe_1d, time_axis, sim_gpe_1d)


@info "computing NPSE" 
if isfile(join([save_path, "npse_heat.jld2"])) && use_precomputed
    @info "\t using precomputed solution npse_heat.jld2" 
    JLD2.@load join([save_path, "npse_heat.jld2"]) npse
else
    sol = runsim(sim_npse; info=false)
    npse = sol.u
    JLD2.@save join([save_path, "npse_heat.jld2"]) npse
end
plot_axial_heatmap(npse, time_axis, sim_npse)


@info "computing NPSE_plus" 
if isfile(join([save_path, "npse_plus_heat.jld2"])) && use_precomputed
    @info "\t using precomputed solution npse_plus_heat.jld2" 
    JLD2.@load join([save_path, "npse_plus_heat.jld2"]) npse_plus
else
    sol = runsim(sim_npse_plus; info=false)
    npse_plus = sol.u
    JLD2.@save join([save_path, "npse_plus_heat.jld2"]) npse_plus
end
plot_axial_heatmap(npse_plus,time_axis, sim_npse_plus, 1)

@info "Done!"

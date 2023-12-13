using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
using DiffEqCallbacks
import CondensateDynamics.V
import FFTW
import JLD2
using Interpolations

gr()
GR.usecolorscheme(1)
include("plot_axial_evolution.jl")
save_path = "results/"

gamma_param = 0.4
use_precomputed = false
maxiters_1d = 100000
N_axial_steps = 512
# =========================================================
## 1D-GPE 

L = (40.0,)
N = (N_axial_steps,)
sim_gpe_1d = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
initial_state = zeros(N[1])

@unpack_Sim sim_gpe_1d
iswitch = -im
equation = GPE_1D
manual = true
solver = SplitStep

# interaction parameter
g_param = gamma_param
maxiters = maxiters_1d
g = - 2 * g_param

n = 100
as = g_param / n
abstol = 1e-6
dt = 0.022
x = X[1]
k = K[1]
dV= volume_element(L, N)
flags = FFTW.EXHAUSTIVE
width = 7
tf = Inf
# SPR condensate bright soliton t in units of omega_perp^-1
analytical_gs = zeros(N)
@. analytical_gs = sqrt(g_param/2) * 2/(exp(g_param*x) + exp(-x*g_param))
psi_0 .= exp.(-(x/1e2).^2)
psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
initial_state .= psi_0
kspace!(psi_0, sim_gpe_1d)
@pack_Sim! sim_gpe_1d

# =========================================================
## NPSE (unable to copy)
L = (40.0,)
N = (N_axial_steps,)
sim_npse = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
initial_state = zeros(N[1])

@unpack_Sim sim_npse
iswitch = -im
equation = NPSE
manual = true
solver = SplitStep

# interaction parameter
g_param = gamma_param
maxiters = maxiters_1d
g = - 2 * g_param

n = 100
as = g_param / n
abstol = 1e-6
dt = 0.022
x = X[1]
k = K[1]
dV= volume_element(L, N)
flags = FFTW.EXHAUSTIVE
width = 7
tf = Inf
psi_0 .= exp.(-(x/1e2).^2)
psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
initial_state .= psi_0
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
initial_state = zeros(N[1])

@unpack_Sim sim_npse_plus
iswitch = -im
equation = NPSE_plus
manual = true
solver = SplitStep

# interaction parameter
g_param = gamma_param
maxiters = maxiters_1d
g = - 2 * g_param

n = 100
as = g_param / n
abstol = 1e-6
dt = 0.022
x = X[1]
k = K[1]
dV= volume_element(L, N)
flags = FFTW.EXHAUSTIVE
width = 7
tf = Inf
psi_0 .= exp.(-(x/1e2).^2)
psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
initial_state .= psi_0
kspace!(psi_0, sim_gpe_1d)
if g_param > 2/3
    @warn "we should expect NPSE collapse"
end
sigma2 = init_sigma2(g)
@pack_Sim! sim_npse_plus

# =========================================================
## 3D-GPE 

L = (40.0,20.0,20.0)
N = (N_axial_steps, 64, 64)
sim_gpe_3d = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)
initial_state = zeros(N[1])

@unpack_Sim sim_gpe_3d
iswitch = -im
equation = GPE_3D
manual = true
solver = SplitStep

g = - g_param * (4*pi)

abstol = 1e-6
maxiters = 100000
dt = 0.5

x0 = 0.0
vv = 0.0

x = Array(X[1])
y = Array(X[2])
z = Array(X[3])
dV= volume_element(L, N)

flags = FFTW.EXHAUSTIVE

tf = Inf
tmp = [exp(-((x-x0)^2+y^2+z^2)/1e6) * exp(-im*x*vv) for x in x, y in y, z in z]
psi_0 = CuArray(tmp)
psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))

kspace!(psi_0, sim_gpe_3d)
tmp = [1/2*(y^2 + z^2) for x in x, y in y, z in z]
V0 = CuArray(tmp)

@pack_Sim! sim_gpe_3d

# @info sim_gpe_3d.g /4/pi
# @info sim_gpe_1d.g /2


# =========================================================
Plots.CURRENT_PLOT.nullableplot = nothing
p = plot_final_density([analytical_gs], sim_gpe_1d; label="analytical", color=:green, doifft=false, ls=:dash)


@info "computing GPE_1D" 
if isfile(join([save_path, "gpe_1d.jld2"])) && use_precomputed
    @info "\t using precomputed solution gpe_1d.jld2" 
    JLD2.@load join([save_path, "gpe_1d.jld2"]) gpe_1d
else
    sol = runsim(sim_gpe_1d; info=false)
    gpe_1d = sol.u
    # JLD2.@save join([save_path, "gpe_1d.jld2"]) gpe_1d
end
plot_final_density!(p, [gpe_1d], sim_gpe_1d; label="GPE_1D", color=:blue, ls=:dash)


@info "computing NPSE" 
if isfile(join([save_path, "npse.jld2"])) && use_precomputed
    @info "\t using precomputed solution npse.jld2" 
    JLD2.@load join([save_path, "npse.jld2"]) npse
else
    sol = runsim(sim_npse; info=true)
    npse = sol.u
    # JLD2.@save join([save_path, "npse.jld2"]) npse
end
plot_final_density!(p, [npse], sim_npse; label="NPSE", color=:blue)


# @info "computing NPSE_plus" 
# if isfile(join([save_path, "npse_plus.jld2"])) && use_precomputed
#     @info "\t using precomputed solution npse_plus.jld2" 
#     JLD2.@load join([save_path, "npse_plus.jld2"]) npse_plus
# else
#     sol = runsim(sim_npse_plus; info=false)
#     npse_plus = sol.u
#     JLD2.@save join([save_path, "npse_plus.jld2"]) npse_plus
# end
# plot_final_density!(p, [npse_plus], sim_npse_plus; label="NPSE_der", ls=:dash, color=:blue)


@info "computing GPE_3D" 
if isfile(join([save_path, "gpe_3d.jld2"])) && use_precomputed
    @info "\t using precomputed solution gpe_3d.jld2" 
    JLD2.@load join([save_path, "gpe_3d.jld2"]) gpe_3d
else
    sol = runsim(sim_gpe_3d; info=true)
    @info size(gpe_3d)
    gpe_3d = sol.u
    # JLD2.@save join([save_path, "gpe_3d.jld2"]) gpe_3d
end

# linear interpolation
gpe_3d = sim_gpe_3d.psi_0
x_axis = sim_npse.X[1] |> real
x_axis_3d = sim_gpe_3d.X[1] |> real
dx = sim_gpe_3d.X[1][2]-sim_gpe_3d.X[1][1]
final_axial = Array(sum(abs2.(xspace(gpe_3d, sim_gpe_3d)), dims=(2, 3)))[:, 1, 1] * sim_gpe_3d.dV / dx |> real
x_3d_range = range(-sim_gpe_3d.L[1]/2, sim_gpe_3d.L[1]/2, length(sim_gpe_3d.X[1])) 
solution_3d = LinearInterpolation(x_3d_range, final_axial, extrapolation_bc = Line())
plot!(p, x_axis, solution_3d(x_axis), label="GPE_3D", color=:black, linestyle=:dot) 


q = plot(x_axis_3d, final_axial, label="GPE_3D", color=:black, linestyle=:dot) 
display(q)
display(p)
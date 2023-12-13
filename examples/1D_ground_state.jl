import GR
using CondensateDynamics
using OrdinaryDiffEq
using DiffEqCallbacks
using LSODA
import CondensateDynamics.V
import FFTW
gr()
GR.usecolorscheme(1)

# =========================================================
## 1D-GPE 
L = (40.0,)
N = (512,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
initial_state = zeros(N[1])

@unpack_Sim sim
iswitch = -im
equation = GPE_1D
manual = true
solver = SplitStep
time_steps = 10000
g_param = 0.5865
g = - 2 * g_param # right??

@info "gamma: " g_param
if equation==NPSE && g_param > 2/3 
    @warn "we should expect NPSE collapse"
end

n = 100
as = - g_param / n
#mu_analytical = npse_mu(n, as)
abstol = 1e-6
maxiters = 10000
dt = 0.022

x = X[1]
k = K[1]
dV= volume_element(L, N)

flags = FFTW.EXHAUSTIVE

tf = Inf
# SPR condensate bright soliton t in units of omega_perp^-1
analytical_gs = zeros(N)
@. analytical_gs = sqrt(g_param/2) * 2/(exp(g_param*x) + exp(-x*g_param))

psi_0 .= exp.(-(x/10).^2)
psi_0 = psi_0 / sqrt(ns(psi_0, sim))
initial_state .= psi_0

if equation == NPSE
    sigma2 = init_sigma2(g) 
end
@info "minimum sigma2" (minimum(sigma2.(psi_0)))
kspace!(psi_0, sim)

@pack_Sim! sim

# =========================================================

sol = runsim(sim; info=false)
final = sol.u

@info "chempot of analytical" chempot(analytical_gs, sim)
@info "final chempot" chempotk(final, sim)
#@info "analytical calculation of chempot" mu_analytical

p = plot(real.(k), abs2.(kspace(initial_state, sim)), color=:blue, ls=:dot, lw=3, label="initial")

plot!(p, real.(k), abs2.(final), color=:red, label="final")
plot!(p, real.(k), abs2.(kspace(analytical_gs, sim)), ls=:dot, lw=2, color=:grey, label="analytical")
display(p)

xspace!(final, sim)

middle = Int(round(N[1]/2))
reduced_x = real.(x)[middle-10:middle+10]
p = plot(real.(x), abs2.(initial_state), color=:green, ls=:dash, lw=2, label="initial")
plot!(p, real.(x), abs2.(final), color=:red, label="final")
plot!(p, real.(x), abs2.(analytical_gs), ls=:dot, lw=2, color=:black, label="analytical")
#plot!(p, real.(x), abs.(V0), ls=:dash, lw=2, color=:green, label="potential")

display(p)
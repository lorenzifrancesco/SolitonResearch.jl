using LaTeXStrings, Plots
import GR
using CondensateDynamics
using OrdinaryDiffEq
import CondensateDynamics.V
gr()
GR.usecolorscheme(1)

## Solve the 1D harmonic oscillator
# problem with 1D-GPE 
L = (40.0,)
N = (256,)
sim = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)

@unpack_Sim sim
g = -2.5
g_param = abs(g) / 2

equation = GPE_1D
manual = true
Nt = 100
t = LinRange(ti,tf,Nt)

iswitch = 1
x = X[1]
k = K[1]
dV = volume_element(L, N)
reltol = 1e-4

x0 = L[1]/4
vv = 6.0
bb = 
tf = x0*2/vv
@. psi_0 = sqrt(g_param/2) * 2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param)) * exp(-im*(x-x0)*vv)
psi_0 = psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = Tsit5()
@. V0 = 4* exp(-10*x^2)

@pack_Sim! sim

sol = runsim(sim; info=false)
final = sol.u[end]
xspace!(final, sim)
xspace!(psi_0, sim)
@info "final distribution norm squared: " ns(final, sim)

p = plot(real.(x), abs2.(psi_0), label="initial")
plot!(p, real.(x), abs2.(final), label="final")
display(p)
xspace!(initial_state, sim)

map(x -> xspace!(x, sim), sol.u[:])
u_mat = mapreduce(permutedims, hcat, sol.u)
ht = heatmap(real.(x), sol.t, abs2.(u_mat)')
display(ht)
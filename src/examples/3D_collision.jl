import JLD2
import CondensateDynamics.V
using CUDA
using LaTeXStrings, Plots
using CUDA.CUFFT
import Makie, GLMakie

JULIA_CUDA_SOFT_MEMORY_LIMIT ="95%"
# ================ plotting functions
function dense(phi)
    psi = xspace(phi,sim)
    density = abs2.(psi)
    pmax = maximum(density)
    if pmax == 0
        throw("Maximum density is null")
    end
    return density/pmax
end

function isosurface_animation(sol,Nt, sim;file="3Dcollision.gif", framerate=3)
    saveto=joinpath("media",file)
    scene = Makie.Scene()
    tindex = Makie.Observable(1)
    iter = [Array(xspace(sol[k], sim)) for k in 1:Nt]
    iter = [abs2.(iter[k]) for k in 1:Nt]
    scene = Makie.volume(Makie.@lift(iter[$tindex]/maximum(iter[$tindex])),
                        algorithm =:iso,
                        isovalue=0.2,
                        isorange=0.1,
                        colormap=:greens, 
                        alpha=0.3,
                        transparency=true
    )

    R = 180
    eyeat = Makie.Vec3f0(R,0,0)
    lookat = Makie.Vec3f0(-50,-50,0)
    Makie.record(scene, saveto, 1:Nt; framerate=framerate) do i
        # Makie.update_cam!(scene, eyeat, lookat)
        # Makie.rotate_cam!(scene, 0., -0.4, 0.)
        tindex[] = i
    end
    return
end

function isosurface(sol)
    scene = Makie.Scene()
    tindex = Makie.Observable(1)
    psol = Array(abs2.(xspace(sol, sim)))
    scene = Makie.volume(psol/maximum(psol),
                        algorithm =:iso,
                        colormap=:greens, 
                        alpha=0.3,
                        isovalue=0.2,
                        isorange=0.1,
    )
    display(scene)
    return
end


# =================== simulation settings
L = (40.0,40.0,40.0)
N = (128,128,128)
sim = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)

# =================== physical parameters 
@unpack_Sim sim
g = -0.587 * 2*pi
g_param = abs(g)/(4*pi)
equation = GPE_3D
iswitch = 1
x = Array(X[1])
y = Array(X[2])
z = Array(X[3])
dV= volume_element(L, N)
reltol = 1e-3
tf = 2
Nt = 30
t = LinRange(ti,tf,Nt)
# nfiles = true
maxiters = 2000
x0 = L[1]/4
vv = 5.0

g_param=3
# tmp = [(exp(-(y^2+z^2)/2) * sqrt(g_param/2)*2/(exp(g_param*(x-x0)) + exp(-(x-x0)*g_param))) * exp(-im*x*vv) for x in x, y in y, z in z]
tmp = [exp(-(y^2+z^2+(x-x0)^2)/2) * exp(-im*x*vv) for x in x, y in y, z in z]

psi_0 = CuArray(tmp)

psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
initial_state = psi_0
kspace!(psi_0, sim)
alg = BS3()
#1/2*(x^2+y^2+ 3*z^2)
tmp = [ (1/2*(y^2+ z^2) + 0.0*1/2*x^2 + 400*exp(-10*x^2)) for x in x, y in y, z in z]
V0 = CuArray(tmp)
#V(x,y,z,t) = 1/2 * (x^2+y^2+z^2)
@pack_Sim! sim


# ===================== simulation
@info "Running solver..."
sol = runsim(sim; info=false)
#final = sol[end]
#JLD2.@save "tmp.jld2" sol
@info "Run complete, plotting..."
# =================== plotting and collect 
#JLD2.@load "tmp.jld2" sol
# gives out of mem
#u = [xspace(sol[k], sim) for k in 1:Nt]

# xspace!(final, sim)
# xspace!(psi_0, sim)
# final = Array(sum(abs2.(sol[end]), dims=(2, 3)))
# psi_0 = Array(sum(abs2.(psi_0), dims=(2, 3)))

# @info "Building animation..."
isosurface_animation(sol,length(sol), sim; framerate=5)
# @info "Completed."
#isosurface(sol[15])
# remember to run 
# sol = nothing
# GC.gc(true)
# to free up GPU mem
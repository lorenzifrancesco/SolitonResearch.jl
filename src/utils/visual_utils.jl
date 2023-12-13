
function single_shot_dynamics(sim::Sim{1, Array{Complex{Float64}}})
    elaps = @elapsed sol = runsim(sim; info=true)
    @info "elapsed time: " elaps
    u = sol.u
    t = sol.t
    @info size(u)
    @info size(t)
    plot_axial_heatmap(u, t, sim; show=true)
    return u
end

function show_psi_0(sim::Sim{1, Array{Complex{Float64}}})
    @unpack N, L, X, psi_0, V0 = sim; x = X[1] |> real
    xpsi_0 = xspace(psi_0, sim)
    p = plot(x, abs2.(xpsi_0), label = "density")
    plot!(p, x, abs.(V0), ls=:dot, color=:grey)
    display(p)
end

function show_psi_0(sim::Sim{3, CuArray{Complex{Float64}}})
    @unpack dV, N, L, X, psi_0, V0 = sim; x = X[1] |> real
    dx = x[2]-x[1]
    xpsi_0 = xspace(psi_0, sim)
    axial_density = sum(abs2.(xpsi_0), dims = (2,3))[:, 1, 1] * dV / dx
    p = plot(x, axial_density, label = "density")
    @warn "not plotting potential... "
    # axial_potential = abs2.(V0)[:, Int(L[2]/2), Int(L[3]/2)]
    # plot!(p, x, axial_potential, ls=:dot, color=:grey)
    display(p)
end

function explain_sim(sim::Sim{1, Array{Complex{Float64}}})
    @unpack_Sim sim
    if iswitch == 1 
        print("\n=== DYNAMICAL SIMULATION ===\n")
    else
        print("\n=== STATICAL  SIMULATION ===\n")
    end
    print("Equation type: $(equation)\n")
    if manual
        print("\t => Manual\n")
    else
        print("\t => Automatic\n")
    end
    print("\n_________Discretization info_________\n")
    print("\t Spatial domain: $(L[1])\n")
    print("\t Spatial points: $(N[1])\n")
    if iswitch == 1
        print("\t Time domain: $(ti) to $(tf)\n")
        print("\t (Time saves: $(Nt))\n")
    else
        print("\t Max iters: $(maxiters)\n")
    end
    print("\n_________Physical info_______________\n")
    print("\t g = $(g)\n")
    print("\t mu = $(mu)\n")
    print("\t Potential maximum = $(maximum(abs.(V0))), minimum = $(minimum(abs.(V0)))\n")
    @pack_Sim
    return nothing
end

function show_sigma2(psi, sim)
    @unpack_Sim sim
    x = X[1] |> real
    sigma2 = estimate_sigma2k(psi, sim)
    p = plot(x, sigma2, label = "σ^2")
    display(p) 
    return nothing
end
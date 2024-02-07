
function single_shot_dynamics(sim::Sim{1,Array{Complex{Float64}}})
    elaps = @elapsed sol = runsim(sim; info = true)
    @info "elapsed time: " elaps
    u = sol.u
    t = sol.t
    @info size(u)
    @info size(t)
    px = plot_axial_heatmap(u, t, sim; show = true)
    savefig(px, "media/checks/single_shot.pdf")
    return u
end

function show_psi_0(sim::Sim{1,Array{Complex{Float64}}})
    @unpack N, L, X, psi_0, V0 = sim
    x = X[1] |> real
    xpsi_0 = xspace(psi_0, sim)
    p = plot(x, abs2.(xpsi_0), label = "density")
    plot!(p, x, abs.(V0), ls = :dot, color = :grey)
    display(p)
    savefig(p, "media/checks/tmp_psi0.pdf")
end

function show_psi_0(sim::Sim{3,CuArray{Complex{Float64}}})
    @unpack dV, N, L, X, psi_0, V0 = sim
    x = X[1] |> real
    dx = x[2] - x[1]
    xpsi_0 = xspace(psi_0, sim)
    axial_density = sum(abs2.(xpsi_0), dims = (2, 3))[:, 1, 1] * dV / dx
    p = plot(x, axial_density, label = "density")
    @warn "not plotting potential... "
    # axial_potential = abs2.(V0)[:, Int(L[2]/2), Int(L[3]/2)]
    # plot!(p, x, axial_potential, ls=:dot, color=:grey)
    display(p)
end

function show_sigma2(psi, sim)
    @unpack_Sim sim
    x = X[1] |> real
    sigma2 = estimate_sigma2k(psi, sim)
    p = plot(x, sigma2, label = "σ^2")
    display(p)
    return nothing
end

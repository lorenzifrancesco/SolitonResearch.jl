using ColorSchemes

function collide(vv = 1.0, bb = 0.0)
    gamma = 0.65
    # add nosaves option
    sd = load_parameters_alt(gamma_param = gamma; eqs = ["G1", "N", "Np"])
    prepare_for_collision!(sd, gamma; use_precomputed_gs = true)
    meas = []
    sim = sd["Np"]
    imprint_vel_set_bar!(sim; vv = vv, bb = bb)
    nn = 2
    pal = palette([:red, :blue], nn)

    mask_refl = map(xx -> xx > 0, sim.X[1] |> real)
    mask_tran = map(xx -> xx < 0, sim.X[1] |> real)
    dt = sim.dt
    if true
        # @unpack_Sim sim
        # dt = 0.001
        # @pack_Sim! sim
        @info "____________________________"
        @info "Computing the transmission"
        @time sol = runsim(sim; info = true)
        @info "Done"
        @info "____________________________"
        final = xspace(sol.u[end], sim)
        pp = plot_final_density(sol.u, sim, label = "final"; show = true, enforce=false)
        ht = plot_axial_heatmap(sol.u, sim.t, sim; show = true, title = "dt=$dt")
        # heatmap(abs2.(sol.u))
        display(ht)
        savefig(ht, "media/checks/tmp_heatmap.pdf")
        savefig(pp, "media/checks/tmp_finals.pdf")
        # readline()
        tran = ns(final, sim, mask_tran)
        refl = ns(final, sim, mask_refl)
        # @assert isapprox(tran + refl, 1.0, atol = 1e-3)
        @info "==> Transmission" tran
        push!(meas, tran)
    else
        dt_range = LinRange(0.01, 1, nn)
        for i = 1:nn
            @unpack_Sim sim
            dt = dt_range[i]
            time_steps = Int(ceil((tf - ti) / dt))
            @warn "dt = " dt
            @warn "tf = " tf
            @pack_Sim! sim
            @time sol = runsim(sim; info = true)
            # plot_axial_heatmap(sol.u, sim.t, sim; show=true, title="dt=$dt")
            # heatmap(abs2.(sol.u))
            # display(ht)
            # readline()

            final = xspace(sol.u[end], sim)
            tran = ns(final, sim, mask_tran)
            refl = ns(final, sim, mask_refl)
            @info "T" tran
            @assert isapprox(tran + refl, 1.0, atol = 1e-3)
            push!(meas, tran)
            GC.gc()
        end
        p = plot(dt_range, meas, color = pal)
        display(p)
    end
    # end
    return meas
end

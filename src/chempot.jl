function compare_chempot(; use_precomputed=true, take_advantage=true)
  pyplot(size=(350, 220))
  sl = load_simulation_list(eqs=[NPSE_plus, GPE_1D, GPE_3D, NPSE])
  N_samples = 50
  gamma_range = LinRange(0.1, 1.0, N_samples) 
  lperp_a = 2e4
  p = plot(xlabel=L"\gamma", ylabel=L"\mu", grid=false)
  q = plot(xlabel=L"\gamma", ylabel=L"E", grid=false, legend=:bottomright, top_margin=5Plots.mm)
  r = plot(xlabel=L"\gamma", ylabel=L"\mathcal{E}", grid=false)
  if isfile("results/mu_db.jld2")
    @info "=> Found dictionary, loading..."
    mud = JLD2.load("results/mu_db.jld2", "mud")
  else
    @info "=> No dictionary found, creating..."
    mud = Dict()
    JLD2.save("results/mu_db.jld2", "mud", mud)
  end
  @info mud

  for sim in sl
    @printf("=========================\n")
    @printf("=== Chempot of [%5s]\n", sim.name)
    @printf("=========================\n")
    sim.iswitch = -im
    mu_vec = zeros(length(gamma_range))
    e_vec = zeros(length(gamma_range))
    e_per_particle_vec = zeros(length(gamma_range))
    k = sim.equation.name
    if haskey(mud, hs(k, 0.666)) && use_precomputed
      mu_vec = mud[hs(k, 0.666)]
    else
      sane = true
      sol = nothing

      for (ig, gamma) in enumerate(gamma_range)
        sim.g = gamma2g(gamma, sim.equation)
        if take_advantage && ig > 1 && sane
          @unpack_Sim sim
          psi_0 .= sol.u[end]
          @pack_Sim! sim
        end
        try
          compound_sol = runsim(sim; info=true)
          sol = compound_sol[1]
          sane = true
          print("-->final time = $(sol.cnt * sim.dt) \n --> Computing chempot...")
          mu_vec[ig] = chempotk(sol.u[end], sim)
        catch e
          if isa(e, NpseCollapse) || isa(e, Gpe3DCollapse)
            if isa(e, NpseCollapse)
              @warn "Collapse detected for gamma = $gamma (NPSE type)"
            else
              @warn "Collapse detected for gamma = $gamma (GPE type)"
            end
            mu_vec[ig] = NaN
            sane = false
            continue
          else
            rethrow(e)
          end
        end
      end
      @warn "pushing"
      push!(mud, hs(k, 0.666) => mu_vec)
      JLD2.save("results/mu_db.jld2", "mud", mud)
    end
    dgamma = gamma_range[2]-gamma_range[1]
    for (ig, gamma) in enumerate(gamma_range)
      if ig == 1
        e_vec[ig] = 0.1
        e_per_particle_vec[ig] = 1.0
      else
        e_vec[ig] = e_vec[ig-1] + (mu_vec[ig]*dgamma)
        e_per_particle_vec[ig] = e_vec[ig] / gamma
      end
    end
    e_vec *= lperp_a
    @info mu_vec
    if sim.equation == NPSE
      plot!(
        p,
        gamma_range,
        clamp!(mu_vec, 0.0, 1.0),
        label=sim.name,
        color=sim.color,
        linestyle=sim.linestyle,
        linewidth=1.5
      )
      plot!(
        q,
        gamma_range,
        e_vec,
        label=sim.name,
        color=sim.color,
        linestyle=sim.linestyle,
        linewidth=1.5
      )
      plot!(
        r,
        gamma_range,
        e_per_particle_vec,
        label=sim.name,
        color=sim.color,
        linestyle=sim.linestyle,
        linewidth=1.5
      )
    else
      plot!(
      p,
      gamma_range,
      clamp!(mu_vec, 0.4, 1.0),
      label=sim.name,
      color=sim.color,
      linestyle=sim.linestyle,
      ylims=(0.4, 1.0)
      )
      plot!(
        q,
        gamma_range,
        e_vec,
        label=sim.name,
        color=sim.color,
        linestyle=sim.linestyle,
        yaxis=(formatter=y->latexstring(round(Int, y * 10^-3))),
      )
      plot!(
        r,
        gamma_range,
        e_per_particle_vec,
        label=sim.name,
        color=sim.color,
        linestyle=sim.linestyle,
      )
    end
  end
  
  annotate!(q, [(0.1, 17800, Plots.text(L"\times10^{3}", 8, :black, :center))])
  savefig(p, "media/chempot_compare.pdf")
  savefig(q, "media/energy_compare.pdf")
  savefig(r, "media/energy_per_particle_compare.pdf")

  plot!(q, xlims=(0.6, 0.8), ylims=(11500, 13900), legend=:topleft)
  annotate!(q, [(0.6047, 14020, Plots.text(L"\times10^{3}", 8, :black, :center))])
  savefig(q, "media/energy_compare_zoom.pdf")

  plot!(r, xlims=(0.6, 0.78), ylims=(2.5, 2.8))
  savefig(r, "media/energy_per_particle_compare_zoom.pdf")
end
function compare_chempot(; use_precomputed=true, take_advantage=true)
  pyplot(size=(350, 220))
  sl = load_simulation_list(eqs=[GPE_1D,NPSE_plus, NPSE, GPE_3D])
  N_samples = 50
  gamma_range = LinRange(0.1, 1.0, N_samples) # TODO
  p = plot(xlabel=L"\gamma", ylabel=L"\mu")
  q = plot(xlabel=L"\gamma", ylabel=L"E")
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
        e_vec[ig] = mu_vec[ig]*dgamma
      else
        e_vec[ig] = e_vec[ig-1]+mu_vec[ig]*dgamma
      end
      end
    plot!(
      p,
      gamma_range,
      mu_vec,
      label=sim.name,
      color=sim.color,
      linestyle=sim.linestyle,
    )
    @info e_vec
    plot!(
      q,
      gamma_range,
      e_vec,
      label=sim.name,
      color=sim.color,
      linestyle=sim.linestyle,
    )
  end
  savefig(p, "media/chempot_compare.pdf")
  savefig(q, "media/energy_compare.pdf")
end
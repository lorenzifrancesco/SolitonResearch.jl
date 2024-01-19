function compare_chempot(; use_precomputed=true, take_advantage=true)
  # pyplot(size=(350, 220))
  sl = load_simulation_list(eqs=[GPE_1D, GPE_3D])
  N_samples = 10
  gamma_range = LinRange(0.0, 1.0, N_samples)
  p = plot(xlabel=L"\gamma", ylabel=L"\mu")

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
    mu_vec = zeros(length(gamma_range))
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
          @time (sol, maxi) = runsim(sim; info=true)
          @warn size(sol.u[end])
          sane = true
          qq = plot()
          plot_final_density!(qq, sol.u, sim; show=true, title=string(gamma))
          display(qq)
          savefig(qq, "media/tmp_gamma_" * string(gamma) * ".pdf")
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
    plot!(
      p,
      gamma_range,
      mu_vec,
      label=sim.equation.name,
      color=:black,
      linestyle=:dashdot,
    )
  end
  savefig(p, "media/chempot_compare.pdf")
end

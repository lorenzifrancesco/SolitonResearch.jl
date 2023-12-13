function compare_chempot(; use_precomputed=true, take_advantage=true)
  pyplot(size=(350, 220))
  @info "Loading parameters..."
  sd = load_parameters_alt(eqs=["G1", "N", "Np","G3"])
  # FIXME
  FFTW.set_num_threads(1)
  gamma_range = LinRange(0.0, 1.0, 40)
  p = plot(xlabel=L"\gamma", ylabel=L"\mu")

  if isfile("results/mu_db.jld2")
    @info "=> Found dictionary, loading..."
    mud = load("results/mu_db.jld2", "mud")
  else
    @info "=> No dictionary found, creating..."
    mud = Dict()
    save("results/mu_db.jld2", "mud", mud)
  end
  @info mud

  for (k, sim) in sd
    @info "==>Running for $k"
    mu_vec = zeros(length(gamma_range))
    
    if haskey(mud, hs(k, 0.666)) && use_precomputed
      @info "--> Loading for $k..."
      mu_vec = mud[hs(k, 0.666)]
    else
      sane = true
      sol = nothing
      @info "--> Computing for $k..."

      for (ig, gamma) in enumerate(gamma_range)
        set_g!(sim, gamma)
        if take_advantage && ig > 1 && sane
          print("\n>>>Taking advantage of previous solution...")
          @unpack_Sim sim
          psi_0 .= sol.u
          @pack_Sim! sim
        end

        print("\n>>>Computing $k for gamma = $gamma")
        try
          sol = runsim(sim; info=true)
          sane = true
          qq = plot_final_density([sol.u], sim; show=true, title=gamma)
          display(qq)
          savefig(qq, "media/tmp.pdf")
          print("-->final time = $(sol.cnt * sim.dt) \n --> Computing chempot...")
          mu_vec[ig] = chempotk(sol.u, sim)
          print("mu = $(mu_vec[ig])\n")
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
      save("results/mu_db.jld2", "mud", mud)
    end
    plot!(p, gamma_range, mu_vec, label=nameof(k), color=colorof(k), linestyle=lineof(k))
  end
  savefig(p, "media/chempot_compare.pdf")
end
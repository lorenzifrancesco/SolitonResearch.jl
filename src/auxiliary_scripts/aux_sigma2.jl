using ColorSchemes

# use the 512 version (long run)
function gs_sigma2(gamma_list=[0.65])
  pyplot(size=(350, 220))
  for gamma in gamma_list
    sd = load_parameters_alt(gamma_param=gamma)
    delete!(sd, "CQ")
    prepare_for_collision!(sd, gamma; use_precomputed_gs=true, info=true)

    p=plot(title="gamma=$(gamma)")
    # if length(sd)>1
    #   pal = palette([:blue, :red], length(sd))
    #   pal = [:blue, :red, :grey, :black]
    #   @assert length(sd) == 4
    # else
    #   pal = [:red]
    # end
    i = 1
    p = plot()
    for (k, v) in sd
      # plot_final_density!(p, [v.psi_0], v; show=true)
      est = estimate_sigma2k(v.psi_0, v)
      if length(v.N) == 3
        est += ones(length(est))* (1 - est[1])
        print("asdvasdvasdv___> ", k )
        plot!(p, real(v.X[1]), circshift(est, Int(ceil(-256/4))), color=colorof(k), linestyle=lineof(k), label=nameof(k))
      else
        plot!(p, real(v.X[1]), circshift(est, Int(ceil(-256/4))), color=colorof(k), linestyle=lineof(k), label=nameof(k))
      end

      if gamma == 0.55
        plot!(p, xlims=(-10, 10), ylims=(0.78, 1.02))
      else
        plot!(p, xlims=(-10, 10), ylims=(0.6, 1.02))
      end
      plot!(grid=false, xlabel=L"x", ylabel=L"\sigma^2")
      # mid = Int(round(v.N[1]/2))
      # mid_tran = Int(round(v.N[2]/2))
      # pp = xspace(v.psi_0, v)
      # tran = real(v.X[2])
      # dx = real(v.X[1][2] - v.X[1][1])
      # dy = real(v.X[2][2] - v.X[2][1])
      # dz = real(v.X[3][2] - v.X[3][1])
      # @warn sum(abs2.(pp[mid, :, mid_tran])/sum(abs2.(pp[mid, :, mid_tran])))
      # q = plot(tran, abs2.(pp[mid, :, mid_tran])/sum(abs2.(pp[mid, :, mid_tran])) / dx)
      # plot!(q, tran, abs2.(pp[mid, mid_tran, :])/sum(abs2.(pp[mid, mid_tran, :])) / dx)
      # plot!(q, tran, exp.(-tran.^2) / sqrt(pi), color=:red, label="Gaussian")
      # plot!(q, tran, exp.(-tran.^2/1.0212) / sqrt(pi*1.0212), ls=:dot, label="Gaussian with correct width")
      # display(q)
      # show_slice(0.0, v.psi_0, v; show=true)
      # show_profile(0.0, v.psi_0, v; show=true)
      i +=1
    end
    plot!(p, legend=:bottomleft)
    # display(p)
    savefig(p, "media/sigma2_"* string(gamma) *".pdf")
  end
end
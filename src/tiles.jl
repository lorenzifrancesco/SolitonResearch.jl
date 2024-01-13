# JULIA_CUDA_SOFT_MEMORY_LIMIT = "95%"

function tiles(;
    extremes = false,
    use_precomputed_gs = true,
    use_precomputed_tiles = false,
    return_maximum = false,
    number_of_tiles = 2,
    equation = "Np",
    messages=false, 
    infos=false,
    plot_finals=false,
    gamma = 0.65
)
    if Threads.nthreads() == 1
        @warn "running in single thread mode!"
    end
    save_path = "results/"

    @info "==============================================="
    @info "\t\tTiling " * string(equation)
    @info @sprintf("\t\tUsing gamma: %.3f", gamma)
    @info "==============================================="
    
    @info "\tLoading parameters..."
    sd = load_parameters(gamma_param = gamma; nosaves = false)
    sd = filter(p -> (first(p) in [equation]), sd)
    @info "\tSetting ground states..."
    @time prepare_for_collision!(sd, gamma; use_precomputed_gs = use_precomputed_gs)

    @info "\tLoading tile library..."
    # create the dictionary
    if isfile(save_path * "tile_dict.jld2")
        @info "\t\tFound library..."
        tile_dict = JLD2.load(save_path * "tile_dict.jld2")
    else
        @info "\t\tNo library! Saving an empty one..."
        tile_dict = Dict()
        JLD2.save(save_path * "tile_dict.jld2", tile_dict)
    end
    sim = sd[equation]
    if haskey(tile_dict, hs(equation, gamma)) && use_precomputed_tiles
        @info "Already found tile for " equation, gamma
    else
        @info "_________________________________________________\n\t Running tiling procedure"
        tile = get_tiles(sim, 
          equation; 
          tiles = number_of_tiles, 
          messages = messages, 
          infos = infos,
          plot_finals = plot_finals
          )
        @info "==== Saving tiles"
        push!(tile_dict, hs(equation, gamma) => tile)
        JLD2.save(save_path * "tile_dict.jld2", tile_dict)
    end
    @info "\tViewing all tiles..."
    plot_tiles()
end


function get_tiles(
    sim::Sim{1,Array{Complex{Float64}}},
    name::String = "noname";
    tiles = 100,
    plot_finals = false,
    messages=true, 
    infos=true
)
    if plot_finals
      @warn "Plotting finals!"
    end
    saveto = "../media/tiles_$(name).pdf"
    max_vel = 1
    max_bar = 1
    #
    vel_list = LinRange(0.1, max_vel, tiles)
    bar_list = LinRange(0, max_bar, tiles)
    tran = Array{Float64,2}(undef, (tiles, tiles))
    refl = Array{Float64,2}(undef, (tiles, tiles))
    warning = zeros((tiles, tiles))

    @info "Filling sim grid..."
    sgrid = Array{Sim,2}(undef, (tiles, tiles))
    archetype = sim
    sgrid[1, 1] = archetype
    @time begin
        for (vx, vv) in enumerate(vel_list)
            for (bx, bb) in enumerate(bar_list)
                sgrid[bx, vx] = imprint_vel_set_bar(archetype; vv = vv, bb = bb)
            end
        end
    end
    # all sims have the same x
    mask_refl = map(xx -> xx > 0, sgrid[1, 1].X[1] |> real)
    mask_tran = map(xx -> xx < 0, sgrid[1, 1].X[1] |> real)

    @info "Running tiling..."
    avg_iteration_time = 0.0
    iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))

    full_time = @elapsed begin
        Threads.@threads for vx in eachindex(vel_list)
        # @showprogress "Computing all the velocities..." for vx in eachindex(vel_list)
            vv = vel_list[vx]
            messages && @info @sprintf("Free memory = %.3f GiB", Sys.free_memory() / 2^30)
            collapse_occured = false
            for (bx, bb) in enumerate(bar_list)
                sim = sgrid[bx, vx]
                @warn sim.dt
                if bx > 2 && isnan(tran[bx-1, vx]) == NaN && isnan(tran[bx-2, vx])
                  messages && @printf("\n Collapse shortcut!")
                  collapse_occured = true
                end
                sol = nothing
                messages && @printf( "\n[tid=%i | bx=%i] (%i/%i)",
                    Threads.threadid(),
                    bx,
                    bx + tiles * (vx - 1),
                    tiles^2
                )
                infos && messages && print("\n")
                if !collapse_occured
                  try
                      avg_iteration_time += @elapsed sol = runsim(sim; info = infos)
                      messages && print("\n --->")
                      # FIXME avoid NPSE+ memory filling problem
                      # @info "GC..."
                      GC.gc()

                      if plot_finals
                          pp = plot_final_density(
                              sol.u,
                              sim;
                              show = false,
                              title = @sprintf("[vx=%i, bx=%i]/%i", vx, bx, tiles)
                          )
                          savefig(pp, "media/checks/final_$(name)_$(vx)_$(bx)_$(tiles).pdf")
                          qq = plot_axial_heatmap(
                              sol.u,
                              sim.t,
                              sim;
                              show = false,
                              title = @sprintf("[vx=%i, bx=%i]/%i", vx, bx, tiles)
                          )
                          savefig(qq, "media/checks/heatmap_$(name)_$(vx)_$(bx)_$(tiles).pdf")
                      end
                  catch err
                      if isa(err, NpseCollapse) || isa(err, Gpe3DCollapse)
                          collapse_occured = true
                      else
                          throw(err)
                      end
                  end
                end
                # catch maxiters hit and set the transmission to zero
                if sim.manual == false
                    if sol.retcode != ReturnCode.Success
                        # @info "Run complete, computing transmission..."
                        @warn "Detected solver failure"
                        tran[bx, vx] = 0.0
                        refl[bx, vx] = 1.0
                        messages && @printf "\t T = %.2f" tran[bx, vx]
                    else
                        # CHANGE : saving the maximum value occured in the iterations
                        final = sol.u[end]
                        # @info "Run complete, computing transmission..."
                        xspace!(final, sim)
                        tran[bx, vx] = ns(final, sim, mask_tran)
                        refl[bx, vx] = ns(final, sim, mask_refl)
                        messages && @printf "\t T = %.2f" tran[bx, vx]
                    end
                else
                    if !collapse_occured
                        final = sol.u[end]
                        # @info "Run complete, computing transmission..."
                        xspace!(final, sim)
                        tran[bx, vx] = ns(final, sim, mask_tran)
                        refl[bx, vx] = ns(final, sim, mask_refl)
                    else
                        messages && print("\tRun complete, detected collapse...")
                        tran[bx, vx] = NaN
                    end
                    messages && @printf "\t T = %.2f" tran[bx, vx]
                end
                if !isapprox(tran[bx, vx] + refl[bx, vx], 1.0, atol = 1e-5)
                    messages &&
                        @printf("\tWARN: [T+R = %.4f]", tran[bx, vx] + refl[bx, vx])
                    warning[bx, vx] = 1.0
                end
            end
        end
    end
    print("\n")
    @info "_________________________________________________"
    @info "Pavement time    = " * @sprintf("%.3f", full_time)
    @info "% time in solver = " * @sprintf("%.3f, %.0f %% of pavement time", avg_iteration_time, avg_iteration_time/full_time*100)
    @info "Single tile time = " * @sprintf("%.3f", avg_iteration_time / tiles^2)
    @info "_________________________________________________"
    print("\n")
    JLD2.@save("tran_$(name).jld2", tran)
    JLD2.@save("refl_$(name).jld2", refl)
    JLD2.@save("warn_$(name).jld2", warning)
    norm_bar = bar_list / max_bar
    norm_vel = vel_list / max_vel
    return tran
end


"""
in the 3D case we do not have sufficient GPU mem, so we go serially
"""
function get_tiles(
    archetype::Sim{3,CuArray{Complex{Float64}}},
    name::String = "noname";
    tiles = 100,
    plot_finals = false,
)
    saveto = "../media/tiles_$(name).pdf"
    max_vel = 1
    max_bar = 1
    #
    vel_list = LinRange(0.1, max_vel, tiles)
    bar_list = LinRange(0, max_bar, tiles)
    tran = Array{Float64,2}(undef, (tiles, tiles))
    refl = Array{Float64,2}(undef, (tiles, tiles))

    @info "Proceeding serially from the archetype..."
    # all sims have the same x
    mask_refl = map(xx -> xx > 0, archetype.X[1] |> real)
    mask_tran = map(xx -> xx < 0, archetype.X[1] |> real)

    @info "Running tiling..."
    avg_iteration_time = 0.0
    iter = Iterators.product(enumerate(vel_list), enumerate(bar_list))
    full_time = @elapsed for ((vx, vv), (bx, bb)) in ProgressBar(iter)
        sim = deepcopy(archetype)
        collapse_occured = false
        imprint_vel_set_bar!(sim; vv = vv, bb = bb)
        @info "Computing tile" (vv, bb)
        sol = nothing
        try
            avg_iteration_time += @elapsed sol = runsim(sim; info = false)
            if plot_finals
                pp = plot_final_density(sol.u, sim; show = false)
                savefig(pp, "media/checks/final_$(name)_$(vv)_$(bb).pdf")
                qq = plot_axial_heatmap(sol.u, sim.t, sim; show = false)
                savefig(qq, "media/checks/heatmap_$(name)_$(vv)_$(bb).pdf")
            end

        catch err
            if isa(err, NpseCollapse) || isa(err, Gpe3DCollapse)
                collapse_occured = true
            else
                throw(err)
            end
        end


        # catch maxiters hit and set the transmission to zero
        if sim.manual == false
            if sol.retcode != ReturnCode.Success
                @info "Run complete, computing transmission..."
                @info "Detected solver failure"
                tran[bx, vx] = 0.0
                refl[bx, vx] = 1.0
                @info "T = " tran[bx, vx]
            else
                if !collapse_occured
                    final = sol.u[end]
                    @info "Run complete, computing transmission..."
                    xspace!(final, sim)
                    tran[bx, vx] = ns(final, sim, mask_tran)
                    refl[bx, vx] = ns(final, sim, mask_refl)
                else
                    @info "Run complete, detected collapse..."
                    tran[bx, vx] = NaN
                end
                @info "T = " tran[bx, vx]
            end
        else
            if !collapse_occured
                final = sol.u[end]
                @info "Run complete, computing transmission..."
                xspace!(final, sim)
                tran[bx, vx] = ns(final, sim, mask_tran)
                refl[bx, vx] = ns(final, sim, mask_refl)
            else
                @info "Run complete, detected collapse..."
                tran[bx, vx] = NaN
                refl[bx, vx] = NaN
            end
            @info "T = " tran[bx, vx]
        end
        if !isapprox(tran[bx, vx] + refl[bx, vx], 1.0, atol = 1e-5)
            @warn "T+R != 1.0"
        end
    end
    @info "Tiling time            = " full_time
    @info "Total time in solver   = " avg_iteration_time
    @info "Average iteration time = " avg_iteration_time / tiles^2

    JLD2.@save("tran_$(name).jld2", tran)
    JLD2.@save("refl_$(name).jld2", refl)
    norm_bar = bar_list / max_bar
    norm_vel = vel_list / max_vel
    return tran
    return tran
end

function plot_tiles()
    # pyplot(size=(300, 220))
    # pyplot(size=(300, 260))
    # backend(:pyplot)
    tile_file = "results/tile_dict.jld2"
    @assert isfile(tile_file)
    td = JLD2.load(tile_file)
    delete!(td, hs("CQ", 0.65))
    ht_list = []
    ct_list = []
    vaxx = []
    baxx = []
    vaxis = nothing
    baxis = nothing
    kk = []
    for (k, v) in td
        push!(kk, k)
        @info "found" ihs(k)
        if ihs(k)[1] == "G3" || ihs(k)[1] == "Np"
            preprocess_tiles_3d!(v)
        end
        baxis = LinRange(0.0, 1.0, size(v)[1])
        vaxis = LinRange(0.1, 1.0, size(v)[1])
        @info v
        push!(vaxx, vaxis)
        push!(baxx, baxis)
        ht2 = heatmap(
            vaxis,
            baxis,
            v,
            clabels = true,
            xlabel = L"v",
            ylabel = L"b",
            colorbar_title = L"T",
            # legend=false,
            aspect_ratio = :equal,
            top_margin = 0 * Plots.mm,
            bottom_margin = 0 * Plots.mm,
            left_margin = 0 * Plots.mm,
            right_margin = 0 * Plots.mm,
        )
        savefig(ht2, "media/tiles_" * string(ihs(k)) * "_ht.pdf")
        push!(ht_list, deepcopy(v))
        v, mask = process_tiles(v)
        ht = contour(vaxis, baxis, v, clabels = true, xlabel = L"v", ylabel = L"b")
        contour!(
            ht,
            vaxis,
            baxis,
            mask,
            levels = [0.0],
            color = :turbo,
            linestyle = :dot,
            linewidth = 1.8,
        )
        push!(ct_list, deepcopy(mask))
        savefig(ht, "media/tiles_" * string(ihs(k)) * "_ct.pdf")
    end

    try
        margins = -4
        heat_first = heatmap(
            vaxx[1],
            baxx[1],
            ht_list[1],
            clabels = true,
            xlabel = L"v",
            ylabel = L"b",
            aspect_ratio = :equal,
            # legend=:none,
            margin = margins * Plots.mm,
        )

        last_idx = length(ht_list)
        heat_last = heatmap(
            vaxx[last_idx],
            baxx[last_idx],
            ht_list[last_idx],
            clabels = true,
            xlabel = L"v",
            aspect_ratio = :equal,
            margin = margins * Plots.mm,
        )

        layout = grid(1, length(ht_list), widths = [0.25, 0.25, 0.25, 0.25])

        ht_comp = plot(
            heat_first,
            [
                heatmap(
                    vaxx[i],
                    baxx[i],
                    ht_list[i],
                    clabels = true,
                    xlabel = L"v",
                    aspect_ratio = :equal,
                    # legend=:none,
                    margin = margins * Plots.mm,
                ) for i in eachindex(ht_list)[2:end-1]
            ]...,
            heat_last,
            layout = layout,
            link = :y,
            leg = false,
            yformatter = _ -> "",
            # size=(900, 250),
            framestyle = :none,
        )

        # display(ht_comp)
        savefig(ht_comp, "media/tiles_ht_comp.pdf")
    catch err
        throw("Not plotting comparison")
    end
end

function preprocess_tiles_3d!(tt)
    for bar = 1:size(tt)[1]
        for vel = 2:size(tt)[2]
            if tt[vel, bar] - tt[vel-1, bar] > 0.02
                tt[vel, bar] = NaN
                for velx = vel:size(tt)[2]
                    tt[velx, bar] = NaN
                end
            end
        end
    end
    return tt
end

function process_tiles(tt)
    mask = ones(size(tt))
    # FIXME the names
    for bar = 1:size(tt)[1]
        for vel = 2:size(tt)[2]
            if abs(tt[vel, bar] - tt[vel-1, bar]) > 0.3
                tt[vel, bar] = NaN
                for velx = 1:vel
                    mask[velx, bar] = 0.0
                end
            end
        end
    end

    flag_prev = true
    flag_curr = true
    posv = 1
    posb = 1
    for vel = 1:size(tt)[2]
        flag_curr = true
        for bar = 1:size(tt)[1]
            if mask[vel, bar] == 0.0
                flag_curr = false
                posv = vel
                posb = bar
            end
        end
        if flag_prev == false && flag_curr == true
            for velx = 1:posv
                mask[velx, posb] = NaN
            end
            return tt, mask
        end
        flag_prev = flag_curr
    end
    return tt, mask
end

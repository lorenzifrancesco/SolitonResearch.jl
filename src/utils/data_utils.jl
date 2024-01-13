
function human_readable_gs(save_path="results/", human_folder="human_readable/")
  if isfile(save_path * "gs_dict.jld2")
    @info "Loading GS library..."
    gs_dict = JLD2.load(save_path * "gs_dict.jld2")
  else
    @info "No GS library found! Quitting."
  end
  display(gs_dict)
  info_file = "infos.txt"
  writedlm(human_folder * info_file, "ciao")
  for (k, gs) in gs_dict
    name = ihs(k)[1]
    gs_file = name * "_complex.dat" ## name, not gamma
    writedlm(human_folder * gs_file, gs)
    gs_file = name * "_abs2.dat" ## name, not gamma
    writedlm(human_folder * gs_file, abs2.(gs))
  end
  nothing
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

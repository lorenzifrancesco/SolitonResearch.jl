
```
 max g allowable for hashing = -5.0, 5.0
```
function hs(eq::String, g::Float64)
  @assert eq in ["G1", "N", "CQ", "Np", "G3"]
  if g <= -5.0
    @warn "Collapse regime selected"
    return string(666666)
  end
  n = 0
  if eq == "G1"
    n += 0
  elseif eq == "N"
    n += 1000
  elseif eq == "Np"
    n += 2000
  elseif eq == "G3"
    n += 3000
  elseif eq == "CQ"
    n += 4000
  else
    throw("Unknown equation")
  end
  n += Int(round(g * 100))
  # print("\nCompute hash: ", n, "\n")
  return string(n)
end


function ihs(s::String)
  n = parse(Int, s)
  if n < 500
    return ("G1", n / 100)
  elseif n < 1500
    return ("N", (n - 1000) / 100)
  elseif n < 2500
    return ("Np", (n - 2000) / 100)
  elseif n < 3500
    return ("G3", (n - 3000) / 100)
  else
    return ("CQ", (n - 4000) / 100)
  end
end

function humanize()
  human_readable_soliton()
  human_readable_tile()
  nothing
end

function human_readable_soliton(
  save_path="results/",
  human_folder="human_readable/",
  file_name="soliton_")
  if isfile(save_path * "gs_dict.jld2")
    @info "Loading GS library..."
    gs_dict = JLD2.load(save_path * "gs_dict.jld2")
  else
    @info "No GS library found! Quitting."
  end
  display(gs_dict)
  info_file = "infos.txt"
  # CSV.write(human_folder * info_file, 4)
  for (k, gs) in gs_dict
    name = ihs(k)[1]
    if name != "G3"
      gs_file = name * "_complex.csv" ## name, not gamma
      CSV.write(human_folder * file_name * gs_file, Tables.table(gs))
      gs_file = name * "_abs2.csv" ## name, not gamma
      CSV.write(human_folder * file_name * gs_file, Tables.table(abs2.(gs)))
    end
  end
  nothing
end

function human_readable_tile(
  save_path="results/",
  human_folder="human_readable/",
  file_name="tiles___")
  if isfile(save_path * "tile_dict.jld2")
    @info "Loading tiles library..."
    tile_dict = JLD2.load(save_path * "tile_dict.jld2")
  else
    @info "No tiles library found! Quitting."
    return
  end
  display(tile_dict)
  for (k, tile) in tile_dict
    name = ihs(k)[1]
    CSV.write(human_folder * file_name * name * ".csv", Tables.table(tile))
  end
  nothing
end

function plot_tiles(matrix=nothing)
  # pyplot(size=(300, 220))
  # pyplot(size=(300, 260))
  # backend(:pyplot)
  tile_file = "results/tile_dict.jld2"
  @assert isfile(tile_file)
  td = JLD2.load(tile_file)
  delete!(td, hs("CQ", 0.65))
  vaxis = nothing
  baxis = nothing
  kk = []
  for (k, v) in td
    push!(kk, k)
    @info "found" ihs(k)
    if ihs(k)[1] == "G3" || ihs(k)[1] == "Np"
      preprocess_tiles_3d!(v)
    end
    (vaxis, baxis) = get_pavement_axes(v)
    plot_pavement(vaxis, baxis, v, k)
  end
  nothing
end

function get_pavement_axes(v)
  baxis = LinRange(0.0, 1.0, size(v)[1])
  vaxis = LinRange(0.1, 1.0, size(v)[1])
  return (vaxis, baxis)
end

function plot_pavement(
  vaxis::LinRange{Float64,Int64},
  baxis::LinRange{Float64,Int64},
  v::Matrix{Float64},
  k::String;
  title=L"T")
  ht2 = heatmap(
    vaxis,
    baxis,
    v,
    clabels=true,
    xlabel=L"v",
    ylabel=L"b",
    colorbar_title=title,
    # legend=false,
    aspect_ratio=:equal,
    top_margin=0 * Plots.mm,
    bottom_margin=0 * Plots.mm,
    left_margin=0 * Plots.mm,
    right_margin=0 * Plots.mm,
  )
  plot!(ht2, vaxis, vaxis .^ 2)
  savefig(ht2, "media/tiles_" * string(ihs(k)) * "_ht.pdf")
  # v, mask = process_tiles(v)
  mask = v
  ht = contour(vaxis, baxis, v, clabels=true, xlabel=L"v", ylabel=L"b")
  contour!(
    ht,
    vaxis,
    baxis,
    mask,
    levels=[0.0],
    color=:turbo,
    linestyle=:dot,
    linewidth=1.8,
  )
  alpha = 1.0
  plot!(ht, vaxis, vaxis .^ 2 * alpha)
  if k == "0000"
    savefig(ht, "media/tiles_test_ct.pdf")
  else
    savefig(ht, "media/tiles_" * string(ihs(k)) * "_ct.pdf")
  end
  nothing
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

function csv2color(file_name=nothing;
  path="results/")
  pyplot()
  if isnothing(file_name)
    file_name = choose_file(path)
    file_name = file_name[1:end-4]
  end
  @warn path*file_name
  vals = Tables.matrix(CSV.read(path * file_name * ".csv", DataFrame))
  (vx, bx) = get_pavement_axes(vals)
  p = heatmap(vx, bx, vals, title=file_name, interpolate=false)
  savefig(p, "media/colorized_csv" * file_name * ".pdf")
  nothing
end

"""
  Choose a file between the ones in a given path.
  Check if the selected number is valid, and if the item is a file
"""
function choose_file(path)
  @printf("Pick a file in <%10s>:\n", path)
  files = readdir(path)
  for (index, file) in enumerate(files)
    println("[$index] $file")
  end
  try
    user_input = parse(Int, readline())
    chosen_file = files[user_input]
    println("You chose: $chosen_file")
    @assert isfile(path*chosen_file)
    return chosen_file
  catch err
    println("Invalid input. Please enter a valid index.")
    throw(err)
  end
end

"""
  Calling git with a shell command, 
  find the commit name and id of the project in pwd
"""
function get_current_commit_data(repo_path=".")
  cmd_cd = `cd $repo_path`
  # run(cmd_cd)
  cmd_git = `git log -1`
  delimiter = "\ngitgitgitgitgitgitgitgitgitgitgitgitgitgitgitgit\n"
  log_string = read(cmd_git, String)
  return log_string
end
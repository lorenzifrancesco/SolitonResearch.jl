function plot_axial_heatmap(
  u,
  time_axis,
  sim::Sim{1,Array{ComplexF64}};
  info=false,
  doifft=true,
  show=false,
  title="",
  path="media/"
)
  @unpack t, X = sim
  x = X[1]
  u = reduce(hcat, u)
  doifft ? u = mapslices(x -> xspace(x, sim), u, dims=(1)) : nothing
  ht = Plots.heatmap(real.(x), t, clamp.(abs2.(u)', -1.0, 2.0), 
  xlabel=L"x", ylabel=L"t", colorbar_title=title, grid=false, minorgrid=false)
  # hidedecorations!(ht, grid=false, minorgrid=false)
  show ? display(ht) : nothing
  savefig(ht, path * "tmp.pdf")
  return ht
end

function plot_axial_heatmap(
  u,
  time_axis,
  sim::Sim{3,CuArray{ComplexF64}};
  axis=1,
  info=false,
  doifft=true,
  show=false,
  title="",)
  print(size(u))
  @unpack t, X = sim
  x = Array(X[axis])
  @assert axis == 1
  ax_list = (1, 2, 3)
  ax_list = filter(x -> x != axis, ax_list)

  doifft ? ux = [xspace(x, sim) for x in u] : nothing
  u_axial = [sum(abs2.(x), dims=ax_list)[:, 1, 1] for x in ux]
  u_axial = Array(reduce(hcat, u_axial))
  ht = Plots.heatmap(real.(x), time_axis, u_axial',  xlabel=L"x", ylabel=L"t",colorbar_title=title)
  display(ht)
  return ht
end


function plot_final_density!(
  p,
  u,
  sim::Sim{1,Array{ComplexF64}};
  info=false,
  doifft=true,
  label="initial",
  lw=1,
  ls=:solid,
  color=:black,
  show=false,
  title=""
  )
  @unpack t, X = sim
  x = Array(X[1])
  tmp = u[end]
  doifft ? final = xspace(tmp, sim) : final = tmp
  @assert isapprox(ns(final, sim), 1.0, atol=1e-3)
  plot!(
    p,
    real.(x),
    abs2.(final),
    label=label,
    linewidth=lw,
    linestyle=ls,
    color=color,
    title=title
  )
  show ? display(p) : nothing
  return p
end

function plot_final_density!(
  p,
  u,
  sim::Sim{3,CuArray{ComplexF64}};
  axis=1,
  info=false,
  doifft=true,
  label="initial",
  show=false,
  lw=1,
  ls=:solid,
  color=:red,
  title="")
  # @error "broken"
  @unpack t, X = sim
  x = Array(X[axis])
  ax_list = (1, 2, 3)
  ax_list = filter(x -> x != axis, ax_list)
  final = u[end]
  doifft ? final = xspace(final, sim) : nothing
  @assert isapprox(ns(final, sim), 1.0, atol=1e-3)
  dydz = (X[2][2]-X[2][1])*(X[3][2]-X[3][1]) |> real
  final_axial = Array(sum(abs2.(final)*dydz, dims=ax_list))[:, 1, 1]
  Plots.plot!(p,
    real.(x),
    final_axial,
    label=label,
    title=title,
    linewidth=lw,
    linestyle=ls,)
  show ? display(p) : nothing
  return p
end

function animation_final_density(u,sim::Sim{1, Array{ComplexF64}};file="1D_evolution.gif",framerate=30,info=false, doifft=true)
    @unpack t, X, Nt = sim; x = Array(X[1]) |> real
    # override until next solution 
    Nt = length(u)
    saveto=joinpath("media",file)
    tindex = Makie.Observable(1)
    iter = u
    doifft ? iter = [Array(xspace(u[k], sim)) for k in 1:Nt] : nothing
    #iter = [ for k in 1:Nt]
    iter = [abs2.(iter[k]) for k in 1:Nt]
    fig, ax= Makie.lines(x, Makie.@lift(iter[$tindex]))

    Makie.limits!(ax, x[1], x[end], 0, 1)
    Makie.record(fig, saveto, 1:Nt; framerate=framerate) do i
        tindex[] = i
    end
    return fig
end
# ================ plotting functions
function dense(phi)
  psi = xspace(phi, sim)
  density = abs2.(psi)
  pmax = maximum(density)
  if pmax == 0
    throw("Maximum density is null")
  end
  return density / pmax
end

function isosurface_animation(u, Nt, sim; file="3Devolution.gif",framerate=3)
    saveto=joinpath("media",file)
    scene = Makie.Scene()
    tindex = Makie.Observable(1)
    iter = [Array(xspace(u[k], sim)) for k in 1:Nt]
    iter = [abs2.(iter[k]) for k in 1:Nt]
    print(size(iter))
    fig = Makie.volume(Makie.@lift(iter[$tindex]/maximum(iter[$tindex])),
                        algorithm =:iso,
                        isovalue=0.2,
                        isorange=0.1,
                        transparency=true
    )

    R = 180
    Makie.record(fig, saveto, 1:Nt; framerate=framerate) do i
        tindex[] = i
    end


    return
end

function transverse_animation(u, Nt, sim; file="section_evolution.gif", fps=4)
  saveto=joinpath("media", "2d"*file)
  midpoint = Int(round(size(u[1])[2]/2))
  iter = [Array(xspace(u[k], sim))[:, :, midpoint] for k in 1:Nt]
  iter = [transpose(abs2.(iter[k])) for k in 1:Nt]
  x = sim.X[1] |> real
  y = sim.X[2] |> real
  anim2d = Plots.Animation()
  maxx = maximum(iter[1])
  for i in 2:Nt
    if maximum(iter[i])>maxx 
      maxx = maximum(iter[i])
    end
  end
  situ_of_max = Array{Int}(undef, sim.N[1])
  for i in 1:Nt
    # y = sin.(x .+ 2π * i / 100)
    # Plots.plot(x, y, ylim=(-1, 1), label="", title="Sine Wave Animation")
    # frame(anim)  # Capture the frame
    situ_of_max[i] = Int(argmax(iter[i][midpoint, :]))
    Plots.heatmap(Array(iter[i]), label="", clims=(0.0, maxx))
    frame(anim2d)
  end
  gif(anim2d, saveto, fps=20)
  
  saveto=joinpath("media", "1d_"*file)
  anim1d = Plots.Animation()

  maxx *= 3
  for i in 1:Nt
    # y = sin.(x .+ 2π * i / 100)
    # Plots.plot(x, y, ylim=(-1, 1), label="", title="Sine Wave Animation")
    # frame(anim)  # Capture the frame
    Plots.plot(abs2.(iter[i][:, situ_of_max[i]]), label="", ylims=(0.0, maxx))
    frame(anim1d)
  end
  gif(anim1d, saveto, fps=20)
end

function isosurface(u)
    scene = Makie.Scene()
    tindex = Makie.Observable(1)
    psol = Array(abs2.(xspace(u, sim)))
    scene = Makie.volume(psol/maximum(psol),
                        algorithm =:iso,
                        isovalue=0.1,
                        isorange=0.1,
    )
    display(scene)
    return
end

function show_slice(
  slice_position::Float64,
  psi,
  sim::Sim{3,CuArray{Complex{Float64}}};
  file="slice.png",
  show=false,
)
  pyplot(size=(250, 200))
  @unpack L, X, N = sim
  psi = Array(xspace(psi, sim))
  x = X[1] |> real
  idx = Int(round((slice_position + L[1] / 2) / L[1] * N[1]))
  @info "idx = $idx"
  @info "L[1]= $(L[1])"
  y = X[2] |> real
  z = X[3] |> real
  saveto = joinpath("media", file)
  ht = heatmap(y, z, abs2.(psi[idx, :, :]), xlabel=L"z", ylabel=L"y", colorbar_title=L"|\psi|^2")
  show ? display(ht) : nothing
  return ht
end

function show_profile(
  slice_position::Float64,
  psi,
  sim::Sim{3,CuArray{Complex{Float64}}};
  file="profile.png",
  show=false,
)
  @unpack L, X, N = sim
  psi = Array(xspace(psi, sim))
  y = X[2] |> real
  idy = Int(round((slice_position + L[2] / 2) / L[2] * N[2]))
  @info "idy = $idy"
  @info "L[2]= $(L[2])"
  x = X[1] |> real
  z = X[3] |> real
  saveto = joinpath("media", file)
  ht = heatmap(x, z, abs2.(psi[:, idy, :])', aspectratio=:equal)
  show ? display(ht) : nothing
  return ht
end

function show_axial_density(
  psi,
  sim::Sim{3,CuArray{Complex{Float64}}};
  file="axial_density.png",
  display=false,
)
  @unpack L, X, N = sim
  psi = Array(xspace(psi, sim))
  x = X[1] |> real
  saveto = joinpath("media", file)
  aa = abs2.(psi)
  dx = x[2] - x[1]
  axial_density = sum(aa, dims=(2, 3))[:, 1, 1] * sim.dV / dx
  p = plot(x, axial_density, label="axial density", color=:red)
  display ? display(p) : nothing
  return p
end

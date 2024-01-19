
function prepare_for_collision!(sim, gamma; 
  use_precomputed_gs=true, 
  info=false,
  verb=false)

  # TODO maybe remove the ground state dictionary, 
  # or make it better with SimSpecsfalse
  verb && @info "______________________________"
  save_path = "results/"
  if isfile(save_path * "gs_dict.jld2")
    verb && @info "Loading GS library..."
    gs_dict = JLD2.load(save_path * "gs_dict.jld2")
  else
    verb && @info "No GS library found! Saving an empty one..."
    gs_dict = Dict()
    JLD2.save(save_path * "gs_dict.jld2", gs_dict)
  end

  name = sim.equation.name
  if haskey(gs_dict, hs(name, gamma)) && use_precomputed_gs
    verb && @info @sprintf("Found in library item (%s, %3.2f)", name, gamma)
  else
    verb && @info @sprintf("Computing item (%s, %3.2f)...", name, gamma)
    uu = get_ground_state(sim; info=info)
    push!(gs_dict, hs(name, gamma) => uu)
    JLD2.save(save_path * "gs_dict.jld2", gs_dict)
  end
  uu = JLD2.load(save_path * "gs_dict.jld2", hs(name, gamma))
  # write the initial state into sim
  if length(sim.N) == 1
    @unpack_Sim sim
    iswitch = 1
    x = X[1]
    psi_0 = uu
    @assert isapprox(nsk(psi_0, sim), 1.0, atol=1e-9)
    # time_steps = Int(ceil((tf - ti) / dt))
    @pack_Sim! sim
  else
    @unpack_Sim sim
    iswitch = 1
    x = X[1] |> real
    y = X[2] |> real
    z = X[3] |> real
    psi_0 = CuArray(uu)
    @assert isapprox(nsk(psi_0, sim), 1.0, atol=1e-9)
    time_steps = Int(ceil((tf - ti) / dt))
    @pack_Sim! sim
  end
  verb && @info @sprintf("Done %s", name)
  verb && @info "______________________________"
  nothing
end

```
Inplace function 1D
```
function imprint_vel_set_bar!(
  sim::Sim{1,Array{Complex{Float64}}};
  vv::Float64=0.0,
  bb::Float64=0.0,
  bw::Float64=0.5,
  dt_set::Float64=0.01,
  time_step_limit::Int64=5000,
  save_each::Bool = false
)

  @unpack_Sim sim
  x = X[1] |> real
  @. V0 = bb * exp(-(x / bw)^2 / 2) # central barrier
  @info @sprintf("barrier width         bw=%8.4f", bw)
  @info @sprintf("barrier crossing time t*=%8.4f", bw/vv)
  @info @sprintf("Nt in crossing period Nt=%8.4f", bw/vv/dt_set)
  x0 = L[1] / 4
  if sim.equation == NPSE_plus
    x0 = L[1] / 8
  else
    x0 = L[1] / 4
  end
  @info x0
  shift = Int(x0 / L[1] * N[1])
  if vv == 0.0
    tf = 2.0
  else
    tf = 2 * x0 / vv
  end
  time_steps = Int(floor((tf - ti) / dt_set))
  if time_steps > time_step_limit
    time_steps = time_step_limit
    dt = (tf - ti) / time_steps
    @warn @sprintf("t_steps > %i, clipped dt=%0.4f", time_step_limit, dt)
  else
    dt = dt_set 
  end
  if save_each
    Nt = time_steps
  else
    Nt = 2
  end
  t = LinRange(ti, tf, Nt)
  xspace!(psi_0, sim)
  psi_0 .= circshift(psi_0, shift)
  @. psi_0 = abs(psi_0) * exp(-im * (x) * vv)
  kspace!(psi_0, sim)
  @pack_Sim! sim
  @warn sim.t
  nothing
end

```
Copy function 1D
```
function imprint_vel_set_bar(
  sim::Sim{1,Array{Complex{Float64}}};
  vv::Float64=0.0,
  bb::Float64=0.0,
  bw::Float64=0.5,
  dt_set::Float64=0.01,
  time_step_limit::Int64=5000,
)
  # @warn "Deepcopy is deprecated"
  simc = deepcopy(sim)
  @unpack_Sim simc
  x = X[1] |> real
  @. V0 = bb * exp(-(x / bw)^2 / 2) # central barrier
  if sim.equation == NPSE_plus
    x0 = L[1] / 8
  else
    x0 = L[1] / 4
  end
  shift = Int(x0 / L[1] * N[1])
  if vv == 0.0
    tf = 2.0
  else
    tf = 2 * x0 / vv
  end
  t = LinRange(ti, tf, Nt)
  time_steps = Int(floor((tf - ti) / dt_set))
  if time_steps > time_step_limit
    time_steps = time_step_limit
    dt = (tf - ti) / time_steps
    # @warn @sprintf("t_steps > %i, clipped dt=%0.4f", time_step_limit, dt)
  end
  xspace!(psi_0, simc)
  psi_0 .= circshift(psi_0, shift)
  @. psi_0 = abs(psi_0) * exp(-im * (x) * vv)
  kspace!(psi_0, simc)
  @pack_Sim! simc
  return simc
end


```
Inplace function 3D
```
function imprint_vel_set_bar!(
  sim::Sim{3,CuArray{Complex{Float64}}};
  vv::Float64=0.0,
  bb::Float64=0.0,
  bw::Float64=0.5,
  dt_set::Float64=0.01, # TODO optimize
  time_step_limit::Int64=5000,
)

  @unpack_Sim sim
  x = X[1] |> real
  y = X[2] |> real
  z = X[3] |> real
  x0 = L[1] / 4
  V0 = [1 / 2 * (z^2 + y^2) + bb * exp(-(x / bw)^2 / 2) for x in x, y in y, z in z] # central barrier
  if vv == 0.0
    tf = 2.0
  else
    tf = 2 * x0 / vv
  end
  shift = Int(x0 / L[1] * N[1])
  t = LinRange(ti, tf, Nt)
  time_steps = Int(floor((tf - ti) / dt_set))
  if time_steps > time_step_limit
    time_steps = time_step_limit
    dt = (tf - ti) / time_steps
    @warn @sprintf("t_steps > %i, clipped dt=%0.4f", time_step_limit, dt)
  end
  xspace!(psi_0, sim)
  psi_0 .= circshift(CuArray(psi_0), (shift, 0, 0))
  @. psi_0 = abs(psi_0) * exp(-im * (x) * vv)
  kspace!(psi_0, sim)
  @pack_Sim! sim
  nothing
end

```
Copy function 3D
```
function imprint_vel_set_bar(
  sim::Sim{3,CuArray{Complex{Float64}}};
  vv::Float64=0.0,
  bb::Float64=0.0,
  bw::Float64=0.5,
  dt_set::Float64=0.01, # TODO optimize
  time_step_limit::Int64=5000,
)
  @warn "Deepcopy is deprecated!"
  simc = deepcopy(sim)
  @unpack_Sim simc
  x = X[1] |> real
  y = X[2] |> real
  z = X[3] |> real
  V0 = [1 / 2 * (z^2 + y^2) + bb * exp(-(x / bw)^2 / 2) for x in x, y in y, z in z] # central barrier
  x0 = L[1] / 4
  if vv == 0.0
    tf = 2.0
  else
    tf = 2 * x0 / vv
  end
  shift = Int(x0 / L[1] * N[1])
  t = LinRange(ti, tf, Nt)
  time_steps = Int(floor((tf - ti) / dt_set))
  if time_steps > time_step_limit
    time_steps = time_step_limit
    dt = (tf - ti) / time_steps
    @warn @sprintf("t_steps > %i, clipped dt=%0.4f", time_step_limit, dt)
  end
  xspace!(psi_0, simc)
  psi_0 .= circshift(CuArray(psi_0), (shift, 0, 0))
  @. psi_0 = abs(psi_0) * exp(-im * (x) * vv)
  kspace!(psi_0, simc)
  @pack_Sim! simc
  return simc
end
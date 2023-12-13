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


function prepare_in_ground_state!(sim::Sim{1, Array{Complex{Float64}}})
    # compute the ground state
    # start from a convenient initial state (it doesn't matter by the way)
    @unpack_Sim sim
    x0 = L[1] / 4
    iswitch = -im 
    maxiters = 1e10
    abstol = 1e-6
    initial_width = 5
    @. psi_0 = exp(-((X[1]-x0)/initial_width)^2/2)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)
    x = X[1] |> real
    V0 = zeros(N[1])
    tmp_dt = deepcopy(dt)
    dt = 0.005
    @pack_Sim! sim

    @info "Computing ground state..."
    sol = runsim(sim; info=false)

    @info "Assigning GS as dynamical sim initial state..."
    xspace!(sol.u, sim)
    display(sol.u)

    # pack everything back up, imprint the correct velocity (suppose x0 stays constant)
    @unpack_Sim sim
    iswitch = 1
    x = X[1]
    dt = tmp_dt
    @. psi_0 = sqrt(abs2(sol.u))
    kspace!(psi_0, sim)
    @assert isapprox(nsk(psi_0, sim), 1.0, atol=1e-9)
    @pack_Sim! sim
    return sim
end

function prepare_in_ground_state!(sim::Sim{3, CuArray{Complex{Float64}}})
    # compute the ground state
    # start from a convenient initial state (it doesn't matter by the way)
    @unpack_Sim sim
    x = X[1] |> real
    y = X[2] |> real
    z = X[3] |> real
    x0 = L[1] / 4
    iswitch = -im 
    maxiters = 1e10
    abstol = 1e-8
    initial_width = 5
    tmp_dt = deepcopy(dt)
    dt = 0.005
    tmp = [exp(-(((x-x0)/initial_width)^2 /2 + (y^2 + z^2)/2)) for x in x, y in y, z in z]
    psi_0 = CuArray(tmp)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim))
    V0 *= 0.0
    psi_0 = kspace(psi_0, sim) # FIXME strangest behaviour ???!
    @pack_Sim! sim

    @info "Computing ground state..."
    sol = runsim(sim; info=true)

    @info "Assigning GS as dynamical sim initial state..."
    xspace!(sol.u, sim)
    print("size of sol.u: ", size(sol.u))

    # pack everything back up, imprint the correct velocity (suppose x0 stays constant)
    @unpack_Sim sim
    iswitch = 1
    x = X[1] |> real
    y = X[2] |> real
    z = X[3] |> real
    dt = tmp_dt
    @. psi_0 = sqrt(abs2.(sol.u))
    psi_0 = psi_0 / sqrt(ns(psi_0, sim))
    kspace!(psi_0, sim)
    @assert isapprox(nsk(psi_0, sim), 1.0, atol=1e-9)
    @pack_Sim! sim
    return sim
end

function prepare_for_collision!(
    sd,
    gamma; 
    use_precomputed_gs=false, 
    info=false
    )
    save_path = "results/"
    if isfile(save_path * "gs_dict.jld2")
        @info "[Loading GS library...]"
        gs_dict = JLD2.load(save_path * "gs_dict.jld2")
        @info "[Done ]"
    else
        @info "No GS library found! Saving an empty one..."
        gs_dict = Dict()
        JLD2.save(save_path * "gs_dict.jld2", gs_dict)
    end
    for (name, sim) in sd
        if haskey(gs_dict, hs(name, gamma)) && use_precomputed_gs
            @info "---> Found in library item " (name, gamma)
        else
            @info "---> Computing item..." (name, gamma)
            uu = get_ground_state(sim; info=info)
            push!(gs_dict, hs(name, gamma) => uu)
            JLD2.save(save_path * "gs_dict.jld2", gs_dict)
        end
        uu = JLD2.load(save_path * "gs_dict.jld2", hs(name, gamma))
        # write the initial state into sim
        @info " ---> Writing ground state into sim..."
        if length(sim.N) == 1
            @unpack_Sim sim
            x0 = L[1] / 4
            shift = Int(x0 / L[1] * N[1])
            iswitch = 1
            x = X[1]
            psi_0 = uu
            xspace!(psi_0, sim)
            psi_0 .= circshift(psi_0, shift)
            kspace!(psi_0, sim)
            @assert isapprox(nsk(psi_0, sim), 1.0, atol=1e-9)
            time_steps  = Int(ceil((tf-ti)/dt))
            @pack_Sim! sim
        else
            @unpack_Sim sim
            x0 = L[1] / 4
            shift = Int(x0 / L[1] * N[1])
            iswitch = 1
            x = X[1] |> real
            y = X[2] |> real
            z = X[3] |> real
            psi_0 = CuArray(uu)
            xspace!(psi_0, sim)
            psi_0 = circshift(CuArray(psi_0), (shift, 0, 0))
            kspace!(psi_0, sim)
            @assert isapprox(nsk(psi_0, sim), 1.0, atol=1e-9)
            time_steps  = Int(ceil((tf-ti)/dt))
            @pack_Sim! sim
        end
    end
    return sd
end

function imprint_vel_set_bar(
    sim::Sim{1, Array{Complex{Float64}}}; 
    vv::Float64=0.0, 
    bb::Float64=0.0,
    bw::Float64=0.5,
    dt::Float64=0.01,
    time_step_limit::Int64=5000)
    simc = deepcopy(sim)
    @unpack_Sim simc
    x = X[1] |> real
    @. V0 = bb * exp(-(x/bw)^2 /2) # central barrier
    x0 = L[1]/4
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    time_steps = Int(floor((tf-ti)/dt))
    if time_steps > time_step_limit
        @warn "time_steps > $time_step_limit, clipping dt"
        time_steps = time_step_limit
        dt = (tf-ti)/time_steps
    end
    xspace!(psi_0, simc)
    @. psi_0 = abs(psi_0) * exp(-im*(x)*vv)
    kspace!(psi_0, simc)
    @pack_Sim! simc
    return simc
end

function imprint_vel_set_bar(
    sim::Sim{3, CuArray{Complex{Float64}}}; 
    vv::Float64=0.0, 
    bb::Float64=0.0, 
    bw::Float64=0.5,
    dt::Float64=0.01, # TODO optimize
    time_step_limit::Int64=5000)

    simc = deepcopy(sim)
    @unpack_Sim simc
    x = X[1] |> real
    y = X[2] |> real
    z = X[3] |> real
    V0 = [1/2*(z^2 + y^2) + bb * exp(-(x/bw)^2 /2) for x in x, y in y, z in z] # central barrier
    x0 = L[1]/4
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    time_steps = Int(floor((tf-ti)/dt))
    if time_steps > time_step_limit
        @warn "time_steps > $time_step_limit, clipping dt"
        time_steps = time_step_limit
        dt = (tf-ti)/time_steps
    end
    xspace!(psi_0, simc)
    @. psi_0 = abs(psi_0) * exp(-im*(x)*vv)
    kspace!(psi_0, simc)
    @pack_Sim! simc
    return simc
end

function imprint_vel_set_bar!(
    sim::Sim{1, Array{Complex{Float64}}}; 
    vv::Float64=0.0, 
    bb::Float64=0.0, 
    bw::Float64=0.5,
    dt::Float64=0.01,
    time_step_limit::Int64=5000)    

    @unpack_Sim sim
    x = X[1] |> real
    @. V0 = bb * exp(-(x/bw)^2 /2) # central barrier
    x0 = L[1]/4
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    time_steps = Int(floor((tf-ti)/dt))
    if time_steps > time_step_limit
        @warn "time_steps > $time_step_limit, clipping dt"
        time_steps = time_step_limit
        dt = (tf-ti)/time_steps
        print("\n Setting dt to $dt")
    end
    xspace!(psi_0, sim)
    @. psi_0 = abs(psi_0) * exp(-im*(x)*vv)
    kspace!(psi_0, sim)
    @pack_Sim! sim
    return sim
end

function imprint_vel_set_bar!(
    sim::Sim{3, CuArray{Complex{Float64}}}; 
    vv::Float64=0.0, 
    bb::Float64=0.0, 
    bw::Float64=0.5,
    dt::Float64=0.01, # TODO optimize
    time_step_limit::Int64=5000)

    @unpack_Sim sim
    x = X[1] |> real
    y = X[2] |> real
    z = X[3] |> real
    x0 = L[1]/4
    V0 = [1/2*(z^2 + y^2) + bb * exp(-(x/bw)^2 /2) for x in x, y in y, z in z] # central barrier
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    time_steps = Int(floor((tf-ti)/dt))
    if time_steps > time_step_limit
        @warn "time_steps > $time_step_limit, clipping dt"
        time_steps = time_step_limit
        dt = (tf-ti)/time_steps
    end
    xspace!(psi_0, sim)
    @. psi_0 = abs(psi_0) * exp(-im*(x)*vv)
    kspace!(psi_0, sim)
    @pack_Sim! sim
    return sim
end

function set_g!(sim::Sim{1, Array{Complex{Float64}}}, gamma_param::Float64=0.4)
    @unpack_Sim sim
    g = -2 * gamma_param
    @pack_Sim! sim
    return
end

function set_g!(sim::Sim{3, CuArray{Complex{Float64}}}, gamma_param::Float64=0.4)
    @unpack_Sim sim
    g = - (4 * pi) * gamma_param
    @pack_Sim! sim
    return
end


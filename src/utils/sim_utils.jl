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


function prepare_for_collision!(sd, gamma; use_precomputed_gs = false, info = false)
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
            uu = get_ground_state(sim; info = info)
            push!(gs_dict, hs(name, gamma) => uu)
            JLD2.save(save_path * "gs_dict.jld2", gs_dict)
        end
        uu = JLD2.load(save_path * "gs_dict.jld2", hs(name, gamma))
        # write the initial state into sim
        @info " ---> Writing ground state into sim..."
        if length(sim.N) == 1
            @unpack_Sim sim
            if sim.equation == "Np" 
              x0 = L[1] / 8
            else
              x0 = L[1] / 4
            end
            shift = Int(x0 / L[1] * N[1])
            iswitch = 1
            x = X[1]
            psi_0 = uu
            xspace!(psi_0, sim)
            psi_0 .= circshift(psi_0, shift)
            kspace!(psi_0, sim)
            @assert isapprox(nsk(psi_0, sim), 1.0, atol = 1e-9)
            time_steps = Int(ceil((tf - ti) / dt))
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
            @assert isapprox(nsk(psi_0, sim), 1.0, atol = 1e-9)
            time_steps = Int(ceil((tf - ti) / dt))
            @pack_Sim! sim
        end
    end
    return sd
end

function imprint_vel_set_bar(
    sim::Sim{1,Array{Complex{Float64}}};
    vv::Float64 = 0.0,
    bb::Float64 = 0.0,
    bw::Float64 = 0.5,
    dt_set::Float64 = 0.001,
    time_step_limit::Int64 = 5000,
)
    simc = deepcopy(sim)
    @unpack_Sim simc
    x = X[1] |> real
    @. V0 = bb * exp(-(x / bw)^2 / 2) # central barrier
    if sim.equation == "Np" 
      x0 = L[1] / 8
    else
      x0 = L[1] / 4
    end
    if vv == 0.0
        tf = 2.0
    else
        tf = 2 * x0 / vv
    end
    t = LinRange(ti, tf, Nt)
    time_steps = Int(floor((tf - ti) / dt_set))
    if time_steps > time_step_limit
        @warn "time_steps > $time_step_limit, clipping dt"
        time_steps = time_step_limit
        dt = (tf - ti) / time_steps
    end
    xspace!(psi_0, simc)
    @. psi_0 = abs(psi_0) * exp(-im * (x) * vv)
    kspace!(psi_0, simc)
    @pack_Sim! simc
    return simc
end

function imprint_vel_set_bar(
    sim::Sim{3,CuArray{Complex{Float64}}};
    vv::Float64 = 0.0,
    bb::Float64 = 0.0,
    bw::Float64 = 0.5,
    dt_set::Float64 = 0.01, # TODO optimize
    time_step_limit::Int64 = 5000,
)

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
    t = LinRange(ti, tf, Nt)
    time_steps = Int(floor((tf - ti) / dt_set))
    if time_steps > time_step_limit
        @warn "time_steps > $time_step_limit, clipping dt"
        time_steps = time_step_limit
        dt = (tf - ti) / time_steps
    end
    xspace!(psi_0, simc)
    @. psi_0 = abs(psi_0) * exp(-im * (x) * vv)
    kspace!(psi_0, simc)
    @pack_Sim! simc
    return simc
end

function imprint_vel_set_bar!(
    sim::Sim{1,Array{Complex{Float64}}};
    vv::Float64 = 0.0,
    bb::Float64 = 0.0,
    bw::Float64 = 0.5,
    dt_set::Float64 = 0.01,
    time_step_limit::Int64 = 5000,
)

    @unpack_Sim sim
    x = X[1] |> real
    @. V0 = bb * exp(-(x / bw)^2 / 2) # central barrier
    x0 = L[1] / 4
    if vv == 0.0
        tf = 2.0
    else
        tf = 2 * x0 / vv
    end
    t = LinRange(ti, tf, Nt)
    time_steps = Int(floor((tf - ti) / dt_set))
    if time_steps > time_step_limit
        @warn "time_steps > $time_step_limit, clipping dt"
        time_steps = time_step_limit
        dt = (tf - ti) / time_steps
        print("\n Setting dt to $dt")
    end
    xspace!(psi_0, sim)
    @. psi_0 = abs(psi_0) * exp(-im * (x) * vv)
    kspace!(psi_0, sim)
    @pack_Sim! sim
    return sim
end

function imprint_vel_set_bar!(
    sim::Sim{3,CuArray{Complex{Float64}}};
    vv::Float64 = 0.0,
    bb::Float64 = 0.0,
    bw::Float64 = 0.5,
    dt_set::Float64 = 0.01, # TODO optimize
    time_step_limit::Int64 = 5000,
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
    t = LinRange(ti, tf, Nt)
    time_steps = Int(floor((tf - ti) / dt_set))
    if time_steps > time_step_limit
        @warn "time_steps > $time_step_limit, clipping dt"
        time_steps = time_step_limit
        dt = (tf - ti) / time_steps
    end
    xspace!(psi_0, sim)
    @. psi_0 = abs(psi_0) * exp(-im * (x) * vv)
    kspace!(psi_0, sim)
    @pack_Sim! sim
    return sim
end

function set_g!(sim::Sim{1,Array{Complex{Float64}}}, gamma_param::Float64 = 0.4)
    @unpack_Sim sim
    g = -2 * gamma_param
    @pack_Sim! sim
    return
end

function set_g!(sim::Sim{3,CuArray{Complex{Float64}}}, gamma_param::Float64 = 0.4)
    @unpack_Sim sim
    g = -(4 * pi) * gamma_param
    @pack_Sim! sim
    return
end

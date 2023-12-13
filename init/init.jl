function load_parameters_gs(; gamma_param::Float64=0.6, eqs=["G1", "N", "Np", "G3"])
    sim_dictionary = Dict()

    maxiters_1d = 1e10
    maxiters_3d = 1e10
    N_axial_steps = 512
    abstol_all = 1e-12
    initial_width = 100 
    dt_all = 0.1
    # =========================================================
    ## 1D-GPE 
    L = (40.0,)
    N = (N_axial_steps,)
    sim_gpe_1d = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
    @unpack_Sim sim_gpe_1d

    iswitch = -im
    equation = GPE_1D
    manual = true
    solver = SplitStep
    # interaction parameter
    maxiters = maxiters_1d
    g = - 2 * gamma_param
    abstol = abstol_all
    dt = dt_all
    x = X[1]
    k = K[1]
    dV= volume_element(L, N)
    flags = FFTW.EXHAUSTIVE
    tf = 1e10
    # SPR condensate bright soliton t in units of omega_perp^-1
    analytical_gs = zeros(N)
    @. analytical_gs = sqrt(gamma_param/2) * 2/(exp(gamma_param*x) + exp(-x*gamma_param))
    @. psi_0 = exp(-x^2/initial_width)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
    kspace!(psi_0, sim_gpe_1d)
    @pack_Sim! sim_gpe_1d

    if "G1" in eqs
        push!(sim_dictionary, "G1" => sim_gpe_1d)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse = deepcopy(sim_gpe_1d)
    initial_state = zeros(N[1])
    @unpack_Sim sim_npse
    equation = NPSE
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse

    if "N" in eqs
        push!(sim_dictionary, "N" => sim_npse)
    end
    # =========================================================
    ## NPSE plus (unable to copy)
    sim_npse_plus = deepcopy(sim_npse)
    @unpack_Sim sim_npse_plus
    equation = NPSE_plus
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse_plus

    if "Np" in eqs
        push!(sim_dictionary, "Np" => sim_npse_plus)
    end
    # =========================================================
    ## 3D-GPE 
    N_axial_steps = 512
    L_axial = 40.0
    L = (L_axial,10.0,10.0)
    N = (N_axial_steps, 64, 64)
    dx = L_axial / N_axial_steps
    sim_gpe_3d = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)
    initial_state = zeros(N[1])
    @unpack_Sim sim_gpe_3d
    iswitch = -im
    equation = GPE_3D
    manual = true
    solver = SplitStep
    g = - gamma_param * (4*pi)
    abstol = abstol_all
    maxiters = maxiters_3d
    dt = dt_all / 10
    x0 = 0.0
    vv = 0.0
    x = Array(X[1])
    y = Array(X[2])
    z = Array(X[3])
    dV= volume_element(L, N)    
    flags = FFTW.EXHAUSTIVE
    tf = 1e10
    tmp = [exp(-((x-x0)^2/initial_width + (y^2 + z^2)/2)) * exp(-im*x*vv) for x in x, y in y, z in z]
    psi_0 = CuArray(tmp)
    psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
    initial_3d = copy(psi_0)
    kspace!(psi_0, sim_gpe_3d)
    tmp = [1/2*(y^2 + z^2) for x in x, y in y, z in z]
    V0 = CuArray(tmp)
    @pack_Sim! sim_gpe_3d
    # @info sim_gpe_3d.g /4/pi
    # @info sim_gpe_1d.g /2

    if "G3" in eqs
        push!(sim_dictionary, "G3" => sim_gpe_3d)
    end
    return sim_dictionary
end

function load_parameters(; vv::Float64 = 0.0, bb::Float64 = 0.0, gamma_param::Float64=0.6, Nsaves::Int64=200, eqs=["G1", "N", "Np", "G3"])
    sim_dictionary = Dict()

    maxiters_1d = 1e10
    maxiters_3d = 1e10
    dt_all = 0.1
    iswitch_all = -im

    max_vel = 1.0
    N_axial_steps = 2048
    abstol_all = 1e-8
    initial_width = 100
    Lx = 40.0
    # =========================================================
    ## 1D-GPE 
    L = (Lx,)
    N = (N_axial_steps,)
    sim_gpe_1d = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
    @unpack_Sim sim_gpe_1d

    iswitch = iswitch_all
    equation = GPE_1D
    manual = true
    solver = SplitStep
    # interaction parameter
    g = - 2 * gamma_param
    n = 100
    abstol = abstol_all
    x = X[1]
    k = K[1]
    x0 = L[1] / 4
    dV= volume_element(L, N)
    flags = FFTW.EXHAUSTIVE
    alg = BS3()

    # will be overwritten
    time_steps = 500
    Nt = Nsaves
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end 
    t = LinRange(ti, tf, Nt)
    dt = (tf-ti)/time_steps

    # specs for GS sim
    maxiters = maxiters_1d
    dt = dt_all

    # SPR condensate bright soliton t in units of omega_perp^-1
    analytical_gs = zeros(N)
    @. analytical_gs = sqrt(gamma_param/2) * 2/(exp(gamma_param*x) + exp(-x*gamma_param)) 
    @. psi_0 = exp(-(x-x0)^2/initial_width) * exp(-im*(x-x0)*vv)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
    kspace!(psi_0, sim_gpe_1d)
    @pack_Sim! sim_gpe_1d

    if "G1" in eqs
        push!(sim_dictionary, "G1" => sim_gpe_1d)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse = deepcopy(sim_gpe_1d)
    initial_state = zeros(N[1])
    @unpack_Sim sim_npse
    equation = NPSE
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse

    if "N" in eqs
        push!(sim_dictionary, "N" => sim_npse)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse_plus = deepcopy(sim_npse)
    @unpack_Sim sim_npse_plus
    equation = NPSE_plus
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse_plus

    if "Np" in eqs
        push!(sim_dictionary, "Np" => sim_npse_plus)
    end
    # =========================================================
    ## 3D-GPE 
    Nx = 512
    L = (Lx,10.0,10.0)
    N = (Nx, 64, 64)
    sim_gpe_3d = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)
    initial_state = zeros(N[1])
    @unpack_Sim sim_gpe_3d
    iswitch = iswitch_all
    equation = GPE_3D
    manual = true
    solver = SplitStep
    g = - gamma_param * (4*pi)
    abstol = abstol_all
    alg = BS3()
    
    x = Array(X[1])
    y = Array(X[2])
    z = Array(X[3])
    dV= volume_element(L, N)    
    flags = FFTW.EXHAUSTIVE
    time_steps = 500
    Nt = Nsaves
    if vv == 0.0
        tf = 2.0
    else
        tf = 2*x0/vv
    end
    t = LinRange(ti, tf, Nt)
    if vv > max_vel/2
        time_steps = 1000
    elseif vv > max_vel/4
        time_steps = 2500
    else
        time_steps = 4000
    end
    dt = (tf-ti)/time_steps

    # specs for GS sim
    maxiters = maxiters_3d
    dt = dt_all

    tmp = [exp(-((x-x0)^2/initial_width + (y^2 + z^2)/2)) * exp(-im*(x-x0)*vv) for x in x, y in y, z in z]
    psi_0 = CuArray(tmp)
    psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
    kspace!(psi_0, sim_gpe_3d)
    tmp = [1/2*(y^2 + z^2) for x in x, y in y, z in z]
    V0 = CuArray(tmp)
    @pack_Sim! sim_gpe_3d
    # @info sim_gpe_3d.g /4/pi
    # @info sim_gpe_1d.g /2

    if "G3" in eqs
        push!(sim_dictionary, "G3" => sim_gpe_3d)
    end
    return sim_dictionary
end

function load_parameters_collapse(
    ; vv::Float64 = 0.0, 
    bb::Float64 = 0.0, 
    gamma_param::Float64=0.6, 
    Nsaves::Int64=200, 
    eqs=["G1", "N", "Np", "G3"],
    initial_width::Float64=5.0,
    )
    sim_dictionary = Dict()

    maxiters_1d = 1e10
    maxiters_3d = 1e10
    dt_all = 0.05
    iswitch_all = -im

    max_vel = 1.0
    N_axial_steps = 128
    abstol_all = 1e-3
    Lx = 60.0
    # =========================================================
    ## 1D-GPE
    L = (Lx,)
    N = (N_axial_steps,)
    sim_gpe_1d = Sim{length(L), Array{Complex{Float64}}}(L=L, N=N)
    @unpack_Sim sim_gpe_1d

    iswitch = iswitch_all
    equation = GPE_1D
    manual = true
    solver = SplitStep
    # interaction parameter
    g = - 2 * gamma_param
    n = 100
    abstol = abstol_all
    x = X[1]
    k = K[1]
    x0 = 0.0 # L[1] / 4
    dV= volume_element(L, N)
    flags = FFTW.EXHAUSTIVE
    alg = BS3()

    # will be overwritten
    time_steps = 500
    Nt = Nsaves
    tf = 2.0
    t = LinRange(ti, tf, Nt)
    dt = (tf-ti)/time_steps

    # specs for GS sim
    maxiters = maxiters_1d
    dt = dt_all

    # SPR condensate bright soliton t in units of omega_perp^-1
    analytical_gs = zeros(N)
    @. analytical_gs = sqrt(gamma_param/2) * 2/(exp(gamma_param*x) + exp(-x*gamma_param)) 
    @. psi_0 = exp(-(x-x0)^2/initial_width) * exp(-im*(x-x0)*vv)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
    kspace!(psi_0, sim_gpe_1d)
    @pack_Sim! sim_gpe_1d
    if "G1" in eqs
        push!(sim_dictionary, "G1" => sim_gpe_1d)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse = deepcopy(sim_gpe_1d)
    initial_state = zeros(N[1])
    @unpack_Sim sim_npse
    equation = NPSE
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse

    if "N" in eqs
        push!(sim_dictionary, "N" => sim_npse)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse_plus = deepcopy(sim_npse)
    @unpack_Sim sim_npse_plus
    equation = NPSE_plus
    # interaction parameter
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse_plus

    if "Np" in eqs
        push!(sim_dictionary, "Np" => sim_npse_plus)
    end
    # =========================================================
    ## 3D-GPE 
    Nx = 512
    L = (Lx,10.0,10.0)
    N = (Nx, 64, 64)
    sim_gpe_3d = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)
    initial_state = zeros(N[1])
    @unpack_Sim sim_gpe_3d
    iswitch = iswitch_all
    equation = GPE_3D
    manual = true
    solver = SplitStep
    g = - gamma_param * (4 * pi)
    abstol = 1e-6
    alg = BS3()
    # we can augment the accuracy


    x = Array(X[1])
    y = Array(X[2])
    z = Array(X[3])
    dV= volume_element(L, N)    
    flags = FFTW.EXHAUSTIVE
    time_steps = 500
    Nt = Nsaves
    tf = 2.0
    t = LinRange(ti, tf, Nt)
    if vv > max_vel/2 # FIXME this should be done later on
        time_steps = 1000
    elseif vv > max_vel/4
        time_steps = 2500
    else
        time_steps = 4000
    end
    dt = (tf-ti)/time_steps

    # specs for GS sim
    maxiters = maxiters_3d
    dt = 0.01

    tmp = [exp(-((x-x0)^2/initial_width + (y^2 + z^2)/2)) * exp(-im*(x-x0)*vv) for x in x, y in y, z in z]
    psi_0 = CuArray(tmp)
    psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
    kspace!(psi_0, sim_gpe_3d)
    tmp = [1/2*(y^2 + z^2) for x in x, y in y, z in z]
    V0 = CuArray(tmp)
    @pack_Sim! sim_gpe_3d
    # @info sim_gpe_3d.g /4/pi
    # @info sim_gpe_1d.g /2

    if "G3" in eqs
        push!(sim_dictionary, "G3" => sim_gpe_3d)
    end
    return sim_dictionary
end


function load_parameters_alt(
    ; vv::Float64 = 0.0,
    bb::Float64 = 0.0, 
    gamma_param::Float64=0.65,
    Nsaves::Int64=200,
    eqs=["G1", "CQ", "N", "Np", "G3"],
    nosaves=false,
    N_axial_1D = 256, ## optimized values
    N_axial_3D = 256,
    N_trans_3D = 40,
    Lx = 40.0,
    Lt = 10.0,
    )

    sim_dictionary = Dict()

    maxiters_1d = 1e10
    maxiters_3d = 1e10
    dt_all = 0.01 # important for the prepare_for_collision function, then overwritten in imprint_vel_set_bar
    iswitch_all = -im

    max_vel = 1.0
    ####
    #### match it to the gs 
    ####
    abstol_all = 1e-3
    # time_steps_all = 200 do not fix it. Use a constant dt

    initial_width = 10
    
    # =========================================================
    ## 1D-GPE 
    L = (Lx,)
    N = (N_axial_1D,)
    sim_gpe_1d = Sim{length(L), Array{Complex{Float64}}}(
      L=L, 
      N=N, 
      )
    @unpack_Sim sim_gpe_1d

    iswitch = iswitch_all
    equation = GPE_1D
    manual = true
    solver = SplitStep
    # interaction parameter
    g = - 2 * gamma_param
    n = 100
    abstol = abstol_all
    x = X[1]
    k = K[1]
    x0 = 0.0 # L[1] / 4
    dV= volume_element(L, N)
    flags = FFTW.EXHAUSTIVE
    alg = BS3()

    # will be overwritten
    if nosaves
      Nt = 2
    else
      Nt = Nsaves
    end
    tf = 2.0
    t = LinRange(ti, tf, Nt)
    dt = dt_all
    time_steps = Int(floor((tf-ti)/dt))
    # specs for GS sim
    maxiters = maxiters_1d

    # SPR condensate bright soliton t in units of omega_perp^-1
    analytical_gs = zeros(N)
    @. analytical_gs = sqrt(gamma_param/2) * 2/(exp(gamma_param*x) + exp(-x*gamma_param)) 
    @. psi_0 = exp(-(x-x0)^2/initial_width) * exp(-im*(x-x0)*vv)
    psi_0 = psi_0 / sqrt(ns(psi_0, sim_gpe_1d))
    kspace!(psi_0, sim_gpe_1d)
    @pack_Sim! sim_gpe_1d

    if "G1" in eqs
        push!(sim_dictionary, "G1" => sim_gpe_1d)
    end

    # # =========================================================
    # ## CQGPE 
    # sim_ccgpe = deepcopy(sim_gpe_1d)
    # @unpack_Sim sim_ccgpe
    # equation = CQGPE
    # @pack_Sim! sim_ccgpe

    # if "CQ" in eqs
    #   push!(sim_dictionary, "CQ" => sim_ccgpe)
    # end

    # =========================================================
    ## NPSE
    sim_npse = deepcopy(sim_gpe_1d)
    initial_state = zeros(N[1])
    @unpack_Sim sim_npse
    equation = NPSE
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse

    if "N" in eqs
        push!(sim_dictionary, "N" => sim_npse)
    end
    # =========================================================
    ## NPSE (unable to copy)
    sim_npse_plus = deepcopy(sim_npse)
    @unpack_Sim sim_npse_plus
    equation = NPSE_plus
    # interaction parameter
    if gamma_param > 2/3
        @warn "we should expect NPSE collapse"
    end
    sigma2 = init_sigma2(g)
    @pack_Sim! sim_npse_plus

    if "Np" in eqs
        push!(sim_dictionary, "Np" => sim_npse_plus)
    end
    # =========================================================
    ## 3D-GPE 
    L = (Lx,Lt,Lt)
    N = (N_axial_3D, N_trans_3D, N_trans_3D)
    sim_gpe_3d = Sim{length(L), CuArray{Complex{Float64}}}(L=L, N=N)
    initial_state = zeros(N[1])
    @unpack_Sim sim_gpe_3d
    iswitch = iswitch_all
    equation = GPE_3D
    manual = true
    solver = SplitStep
    g = - gamma_param * (4 * pi)
    abstol = abstol_all
    alg = BS3()
    # we can augment the accuracy


    x = Array(X[1])
    y = Array(X[2])
    z = Array(X[3])
    dV= volume_element(L, N)    
    flags = FFTW.EXHAUSTIVE
    if nosaves
      Nt = 2
    else
      Nt = 50
    end
    tf = 2.0
    t = LinRange(ti, tf, Nt)
    dt = dt_all
    time_steps = Int(floor((tf-ti)/dt))

    # specs for GS sim
    maxiters = maxiters_3d

    tmp = [exp(-((x-x0)^2/initial_width + (y^2 + z^2)/2)) * exp(-im*(x-x0)*vv) for x in x, y in y, z in z]
    psi_0 = CuArray(tmp)
    psi_0 .= psi_0 / sqrt(sum(abs2.(psi_0) * dV))
    kspace!(psi_0, sim_gpe_3d)
    tmp = [1/2*(y^2 + z^2) for x in x, y in y, z in z]
    V0 = CuArray(tmp)
    @pack_Sim! sim_gpe_3d
    # @info sim_gpe_3d.g /4/pi
    # @info sim_gpe_1d.g /2

    if "G3" in eqs
        push!(sim_dictionary, "G3" => sim_gpe_3d)
    end
    return sim_dictionary
end
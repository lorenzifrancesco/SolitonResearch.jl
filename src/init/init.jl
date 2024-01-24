"""
Load the simulation dictionary.
Input: directory of input files
Output: dictionary with selected 
"""
function load_simulation_list(;
  input_dir="input/",
  eqs=[GPE_3D, GPE_1D, NPSE, NPSE_plus],
  verb=true,
  idx_domain=1)
  sim_dictionary::Array{Sim} = []
  for eq in eqs
    push!(sim_dictionary, load_simulation(input_dir, 
                                          eq, 
                                          verb=verb,
                                          idx_domain=idx_domain))
  end
  #  sort(sim_dictionary, lt=SolitonDynamics.isless)
  return sim_dictionary
end

"""
Load a single simulation from the input specification
Input: input configuration directory, equation selection, row of the configuration file to be loaded
"""
function load_simulation(input_dir, eq::EquationType;
  idx_domain=1,
  idx_nonlin=1,
  idx_precis=1,
  verb=false
)
  domain_df = CSV.read(input_dir * "domain.csv", DataFrame)
  nonlin_df = CSV.read(input_dir * "nonlinearity.csv", DataFrame)
  precis_df = CSV.read(input_dir * "precision.csv", DataFrame)
  verb && @info @sprintf("Loading %s... \n\tDomain      : %s,\n\tNonlinearity: %s,\n\tPrecis      : %s.", eq.name, domain_df.name[idx_domain], nonlin_df.name[idx_nonlin], precis_df.name[idx_precis])
  display(domain_df)
  display(nonlin_df)
  display(precis_df)
  if eq == GPE_3D
    L = (domain_df.L_axial[idx_domain], domain_df.L_radial[idx_domain], domain_df.L_radial[idx_domain])
    N = (domain_df.N_axial_3D[idx_domain], domain_df.N_radial[idx_domain], domain_df.N_radial[idx_domain])
    sim = Sim{length(L),CuArray{Complex{Float64}}}(L=L, N=N)
  else
    L = (domain_df.L_axial[idx_domain],)
    N = (domain_df.N_axial_1D[idx_domain],)
    sim = Sim{length(L),Array{Complex{Float64}}}(L=L, N=N)
  end
  sim.name = eq.name
  sim.iswitch = 1.0
  sim.manual = true
  sim.equation = eq
  sim.solver = SplitStep
  # interaction parameter
  gamma_param = nonlin_df.gamma[idx_nonlin]
  sim.g = gamma2g(gamma_param, sim.equation)
  if eq == GPE_3D
    check = -gamma_param * (4 * pi)
  else
    check = -2 * gamma_param
  end

  @assert check == sim.g
  if eq in [NPSE, NPSE_plus]
    if gamma_param > 2 / 3
      @warn "we should expect NPSE collapse"
    end
    sim.sigma2 = init_sigma2(sim.g)
  end
  sim.dV = volume_element(L, N)
  # sim.flags = FFTW.EXHAUSTIVE
  sim.collapse_threshold = nonlin_df.collapse_threshold[idx_nonlin]
  sim.abstol = precis_df.abstol[idx_precis]
  if eq == GPE_3D
    x = Array(sim.X[1])
    y = Array(sim.X[2])
    z = Array(sim.X[3])
  else
    x = sim.X[1]
  end
  if precis_df.no_saves[idx_precis]
    sim.Nt = 2
  else
    sim.Nt = precis_df.N_saves[idx_precis]
  end
  sim.tf = 2.0
  sim.t = LinRange(0.0, sim.tf, sim.Nt)
  sim.dt = domain_df.dt[idx_domain]
  sim.time_steps = Int(floor(sim.tf / sim.dt))
  # specs for GS sim

  # SPR condensate bright soliton t in units of omega_perp^-1
  initial_width = nonlin_df.initial_width[idx_nonlin]
  if eq == GPE_3D
    tmp = [exp(-((x)^2 / initial_width + (y^2 + z^2) / 2)) for x in x, y in y, z in z]
    sim.psi_0 .= CuArray(tmp)
    sim.psi_0 .= sim.psi_0 / sqrt(sum(abs2.(sim.psi_0) * sim.dV)) # is this normalization needed because of CuArray/Array? 
    sim.maxiters = domain_df.max_iters[idx_domain] / 10 # set smart maxiters in 3D
    tmp .= [1 / 2 * (y^2 + z^2) for x in x, y in y, z in z]
    sim.V0 = CuArray(tmp)
  else
    @. sim.psi_0 = exp(-(x / initial_width)^2)
    sim.psi_0 = sim.psi_0 / sqrt(ns(sim.psi_0, sim))
    sim.maxiters = domain_df.max_iters[idx_domain]
  end
  kspace!(sim.psi_0, sim)
  @assert isapprox(nsk(sim.psi_0, sim), 1.0, rtol=1.0e-9)
  sim.color = get_color(sim.equation)
  sim.linestyle = get_linestyle(sim.equation)
  sim.name = paper_name(sim.equation)
  return sim
end

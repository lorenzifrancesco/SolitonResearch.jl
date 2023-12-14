function efficiency(
  select_eq = "Np",
  vv = 1.0,
  bb = 1.0,
  dt = 0.1
  )

  sd = load_parameters_alt()
  n = sd["N"]
  np = sd["Np"]
  xspace!(n.psi_0, n)
  kspace!(np.psi_0, np)
  xspace!(np.psi_0, np)
  @info "done"
  prepare_for_collision!(sd, 0.65; use_precomputed_gs=true)
  # sim = sd[select_eq]
  # imprint_vel_set_bar!(sim; vv=vv, bb=bb, dt=dt)
  # @info sim
  # @info sim.time_steps
  # @time single_shot_dynamics(sim)
  return
end

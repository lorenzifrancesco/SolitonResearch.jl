using SolitonDynamics, SolitonResearch, CUDA
CUDA.allowscalar()
sim = load_simulation("input/", GPE_3D)
prepare_for_collision!(sim, 0.65; use_precomputed_gs=false)
sim.g /=3
sol = runsim(sim)

u = sol[1].u
transverse_animation(u, sim.Nt, sim)
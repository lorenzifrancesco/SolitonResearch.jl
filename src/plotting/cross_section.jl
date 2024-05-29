using SolitonDynamics, SolitonResearch, CUDA
CUDA.allowscalar()
sim = load_simulation("input/", GPE_3D)
prepare_for_collision(sim, 0.67*3; use_precomputed=false)
sol = runsim(sim)
u = sol[1].u
transverse_animation(u, sim.Nt, sim)
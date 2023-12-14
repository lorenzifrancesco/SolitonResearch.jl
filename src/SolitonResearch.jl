module SolitonResearch 

using ExportAll
using PrecompileTools

using SolitonDynamics, CUDA, FFTW, OrdinaryDiffEq
using Distributed
using Plots
using HDF5
import JLD2
using Interpolations
using OrderedCollections
using Printf
using LaTeXStrings
# import Makie, GLMakie
using ProgressBars, Colors, ColorSchemes

# FFTW.set_num_threads(2)

include("init/_plot_settings.jl")
include("plotting/plot_axial_evolution.jl")
include("plotting/plot_isosurfaces.jl")
include("utils/visual_utils.jl")
include("init/init.jl")
include("utils/sim_utils.jl")
include("solitons.jl")
include("lines.jl")
include("tiles.jl")
include("chempot.jl")
include("efficiency.jl")
include("auxiliary_scripts/aux_collapse.jl")
include("auxiliary_scripts/aux_gs.jl")
include("auxiliary_scripts/aux_collision.jl")
include("auxiliary_scripts/aux_sigma2.jl")

@exportAll()

# @setup_workload begin
#   @compile_workload begin
#     @info "entering compile workload"
#     sd = load_parameters_alt()
#     # maybe too much
#     # prepare_for_collision!(sd, 0.65)
#   end
# end

end
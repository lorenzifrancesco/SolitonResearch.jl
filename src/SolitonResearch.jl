module SolitonResearch

using ExportAll
using PrecompileTools

using SolitonDynamics
import CUDA: CuArray 
using Plots
import JLD2, CSV
using Tables, DataFrames
using Interpolations
using OrderedCollections
using Printf
using LaTeXStrings
using ProgressBars, ProgressMeter, Colors, ColorSchemes 
using Distributed
 
include("init/_plot_settings.jl")
include("plotting/plot_axial_evolution.jl")
include("plotting/plot_isosurfaces.jl")
include("utils/visual_utils.jl")
include("init/init.jl")
include("utils/data_utils.jl")
include("utils/sim_utils.jl")
include("solitons.jl")
include("lines.jl")
include("tiles.jl")
include("chempot.jl")

@exportAll()
export GPE_1D, NPSE, NPSE_plus, GPE_3D

# @setup_workload begin
#   @compile_workload begin
#     @info "entering compile workload"
#     sd = load_parameters()
#     # prepare_for_collision!(sd, 0.65)
#   end
# end

end

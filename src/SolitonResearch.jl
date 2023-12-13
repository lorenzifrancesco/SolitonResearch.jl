module SolitonResearch 

using ExportAll
using PrecompileTools

using SolitonDynamics, CUDA, FFTW, OrdinaryDiffEq
# using Plots
using HDF5
import JLD2
using Interpolations
using OrderedCollections

using LaTeXStrings
import Makie, GLMakie
using ProgressBars, Colors, ColorSchemes

# using PyCall
# const plt = pyimport("matplotlib.pyplot")

# includet("/home/lorenzi/SolitonDynamics.jl/src/CondensateDynamics.jl")
# using Main.CondensateDynamics

# # Set other parameters as needed
# plt.rcParams["figure.figsize"] = [8, 6]
# plt.rcParams["font.size"] = 12
# plt.rcParams["lines.linewidth"] = 2
# Set other parameters as needed

pyplot(size=(350, 220))

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
include("auxiliary_scripts/aux_collapse.jl")
include("auxiliary_scripts/aux_gs.jl")
include("auxiliary_scripts/aux_collision.jl")
include("auxiliary_scripts/aux_sigma2.jl")

@exportAll()

@setup_workload begin
  @compile_workload begin
    @info "entering compile workload"
    sd = load_parameters_alt()
    prepare_for_collision!(sd, 0.65)
  end
end

end
module SolitonResearch 

using ExportAll

using SolitonDynamics
using PyPlot
using Revise
using HDF5
import JLD2
using FFTW, CUDA, OrdinaryDiffEq
using Interpolations
using OrderedCollections
using LoopVectorization, LinearAlgebra

# ENV["PYTHON"]="usr/bin/python"
# using PyPlot
using LaTeXStrings
using Plots
import Makie, GLMakie
using ProgressBars, Colors

pyplot()
using PyCall
const plt = pyimport("matplotlib.pyplot")

# includet("/home/lorenzi/SolitonDynamics.jl/src/CondensateDynamics.jl")
# using Main.CondensateDynamics

# # Set other parameters as needed
# plt.rcParams["figure.figsize"] = [8, 6]
# plt.rcParams["font.size"] = 12
# plt.rcParams["lines.linewidth"] = 2
# Set other parameters as needed

includet("init/_plot_settings.jl")
# pyplot(size=(350, 220))
# if ENV["USER"] == "ubuntu"
#   plotly()
# else
#   pyplot()
# end

includet("plotting/plot_axial_evolution.jl")
includet("plotting/plot_isosurfaces.jl")
includet("utils/visual_utils.jl")
includet("init/init.jl")
includet("utils/sim_utils.jl")

includet("solitons.jl")
includet("lines.jl")
includet("tiles.jl")
includet("chempot.jl")
includet("auxiliary_scripts/aux_collapse.jl")
includet("auxiliary_scripts/aux_gs.jl")
includet("auxiliary_scripts/aux_collision.jl")
includet("auxiliary_scripts/aux_sigma2.jl")
@exportAll
end
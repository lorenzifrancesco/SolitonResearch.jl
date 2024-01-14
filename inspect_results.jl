using SolitonResearch, CSV
using Plots
using Tables, DataFrames
pyplot()
results_dir = "results/"
tran = Tables.matrix(CSV.read(results_dir*"tran.csv", DataFrame))
(vx, bx) = get_pavement_axes(tran)
display(tran)
plot_pavement(vx, bx, tran, "0000")
using CSV, DataFrames, Tables, Printf
using SolitonResearch

# conf = CSV.read("conf.csv", DataFrame)
# @info conf
# @info @sprintf("LAUNCHING %ix%i tiles, equation = %s", conf.N[1], conf.N[1], conf.eq[1])
N = 10
return_maximum = true
@info @sprintf("N=%i, return_maximum=%d", N, return_maximum)
fill_tiles(number_of_tiles=N, eqs=[NPSE_plus], return_maximum=return_maximum)
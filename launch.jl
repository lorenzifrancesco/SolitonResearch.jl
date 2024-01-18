using CSV, DataFrames, Tables, Printf
using SolitonResearch

conf = CSV.read("conf.csv", DataFrame)
@info conf
@info @sprintf("LAUNCHING %ix%i tiles, equation = %s", conf.N[1], conf.N[1], conf.eq[1])
fill_tiles(number_of_tiles=50, eqs=[NPSE_plus], return_maximum=true)

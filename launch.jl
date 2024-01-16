using CSV, DataFrames, Tables, Printf
using SolitonResearch

conf = CSV.read("conf.csv", DataFrame)
@info conf
@info @sprintf("LAUNCHING %ix%i tiles, equation = %s", conf.N[1], conf.N[1], conf.eq[1])
tiles(number_of_tiles=50, equation="G1", messages=true)
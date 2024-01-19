using CSV, DataFrames, Tables, Printf
using SolitonResearch

<<<<<<< HEAD
conf = CSV.read("conf.csv", DataFrame)
@info conf
@info @sprintf("LAUNCHING %ix%i tiles, equation = %s", conf.N[1], conf.N[1], conf.eq[1])
fill_tiles(number_of_tiles=50, eqs=[NPSE_plus], return_maximum=true)
=======
print("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@ BEGIN OF RUN @@@@@@@@@@@@@@@@@@@@@@@@@@\n")
cmd = `date`
print(read(cmd, String))
print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

# TODO load configuration from input
N = 10
return_maximum = true
@printf("N=%i, return_maximum=%d", N, return_maximum)
print("\n=================================================================\n")
print(get_current_commit_data())
print("=================================================================\n")
fill_tiles(number_of_tiles=N, eqs=[NPSE_plus], return_maximum=return_maximum)
print("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF RUN @@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
>>>>>>> 56ab84b491822295de54da76f93cc6b60d67f63b

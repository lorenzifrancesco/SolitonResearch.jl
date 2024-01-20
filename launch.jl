using CSV, DataFrames, Tables, Printf
using SolitonResearch, Plots; gr()

print("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@ BEGIN OF RUN @@@@@@@@@@@@@@@@@@@@@@@@@@\n")
cmd = `date`
print(read(cmd, String))
print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

# TODO load configuration from input
print("\n=================================================================\n")
print(get_current_commit_data())
print("=================================================================\n")
N = 20
return_maximum = true
@printf("==============      N=%i, return_maximum=%d       ================\n", N, return_maximum)
fill_tiles(number_of_tiles=N, eqs=[NPSE_plus], return_maximum=return_maximum)

print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
cmd = `date`
print(read(cmd, String))
print("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF RUN @@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
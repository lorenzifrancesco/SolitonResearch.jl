using Printf
print("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@ BEGIN OF RUN @@@@@@@@@@@@@@@@@@@@@@@@@@\n")
cmd = `date`
print(read(cmd, String))
print("NAME: ")
@printf("[%30s]", readline())

using CSV, DataFrames, Tables, Printf
using SolitonResearch, Plots; gr()

print("\n@@@ packages are loaded @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

# TODO load configuration from input
print("\n=================================================================\n")
print(get_current_commit_data())
print("\n=================================================================\n")
N = 50
return_maximum = true
@printf("==============      N=%i, return_maximum=%d       ================\n", N, return_maximum)
fill_tiles(number_of_tiles=N, eqs=[GPE_1D, NPSE], return_maximum=return_maximum)

print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
cmd = `date`
print(read(cmd, String))
print("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF RUN @@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
#! /bin/bash
julia --threads=1 launch.jl > >(tee -a stdout_20.log) 2> >(tee -a stderr_20.log >&2)



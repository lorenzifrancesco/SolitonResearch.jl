#! /bin/bash
julia --threads=1 launch.jl > >(tee -a stdout_50.log) 2> >(tee -a stderr_50.log >&2)
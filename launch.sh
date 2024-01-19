#! /bin/bash
julia --threads=1 launch.jl > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)


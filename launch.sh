#! /bin/bash
date | tee -a stdout.log
date >> stderr.log
julia --threads=12 launch.jl > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)


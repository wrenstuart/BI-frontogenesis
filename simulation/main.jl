# This is all written for Oceananigans v0.90.6
# Note: 512×512×64 sims take up 6247MiB of memory

const n_d = 20      # for n_d×n_d array of drifters
include("sim.jl")
include("inputs/" * ARGS[1] * ".jl")

run_sim(sim_params(), ARGS[1])
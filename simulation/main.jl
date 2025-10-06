# Define sets of Richardson numbers (Ri) and isopycnal slopes (α) over which to iterate
# These numbers are defined by their base-10 logarithms

const n_d = 20      # for n_d×n_d array of drifters
include("sim.jl")
include("inputs/" * ARGS[1] * ".jl")

run_sim(sim_params(), ARGS[1])
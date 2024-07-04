# Define sets of Richardson numbers (Ri) and isopycnal slopes (α) over which to iterate
# These numbers are defined by their base-10 logarithms

include("sim-ptcl-test-2.jl")
include(ARGS[1] * "_input.jl")

run_sim(sim_params())
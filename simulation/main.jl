# This is all written for Oceananigans v0.90.6
# Note: 512×512×64 sims take up 6247MiB of memory

const n_d = 20      # for n_d×n_d array of drifters
label::String = ARGS[1]
include("sim.jl")
include("inputs/" * label * ".jl")

@info label
dir = "raw_data/" * label * "/"
if isdir(dir)
    if length(readdir(dir)) > 0
        throw("Output directory for label " * label * " already used")
    end
    @info "Outputting to empty directory: " * dir
else
    mkdir(dir)
    @info "Created directory for simulation with label " * label
end
run_sim(sim_params(), label)
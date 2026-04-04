# This is all written for Oceananigans v0.90.6
# Note: 512×512×64 sims take up 6247MiB of memory

const n_d = 20      # for n_d×n_d array of drifters
label::String = ARGS[1]
include("sim.jl")
include("inputs/" * label * ".jl")
include("../io.jl")

function containscheckpoint(dir::String)
    for filename in readdir(dir)
        if occursin("checkpoint", filename)
            return true
        end
    end
    return false
end

dir = "raw_data/" * label * "/"
if isdir(dir)
    if containscheckpoint(dir)
        run_sim(sim_params(), label; pickup = true)
    end
    if length(readdir(dir)) > 0
        throw("Output directory for label " * label * " already used")
    end
    print("Outputting to empty directory: " * dir, label)
else
    mkdir(dir)
    print("Created directory for simulation with label " * label, label)
end
touch(dir * "log.txt")
doubleoutput("Log file for simulation with label " * label, label)
run_sim(sim_params(), label)
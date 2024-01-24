To use, run main.jl

First, must create a small file simulation/simulation_parameters.jl, containing one function:

function sim_params()
    par_with_GPU = false    # Whether you want to parallelise on a CPU or GPU
    res = (128, 128, 16)    # The resolution of the simulation
    return (GPU = par_with_GPU, res = res)
end

This is not tracked by git because it will vary between machines.
Plan to eventually move log(Ri) and log(s) into this file
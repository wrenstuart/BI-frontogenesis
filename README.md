To use, run main.jl

First, must create a small file simulation/[name]_input.jl, containing one function, sim_params. Here is an example:

function sim_params()
    label = "test"
    Ri = 1e3
    s = 1e4
    ν_h = 1e-4
    ν_v = 1e-3
    par_with_GPU = false
    res = (128, 128, 16)
    return (GPU = par_with_GPU, res = res, Ri = Ri, s = s, ν_h = ν_h, ν_v = ν_v, label = label)
end

This is not tracked by git.

To run the simulation, use julia to run
    main.jl [name]
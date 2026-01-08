# Now outputting halo points, which we would like to check
function sim_params() 
    Ri = 1
    s = 1e4
    ν_v = 1e-3
    ν_h = 4e+1
    GPU = true
    res = (10, 10, 10)
    advection_scheme = Centered
    horizontal_hyperviscosity = false
    short_duration = true
    diffusive_cfl = 0.1
    return (; GPU, res, Ri, s, ν_h, ν_v, advection_scheme, horizontal_hyperviscosity, short_duration, diffusive_cfl)
end
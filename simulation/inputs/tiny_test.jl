function sim_params() 
    Ri = 1
    s = 1e4
    ν_v = 1e-2
    ν_h = 4e+2
    GPU = false
    res = (32, 32, 8)
    advection_scheme = CenteredSecondOrder
    horizontal_hyperviscosity = false
    short_duration = true
    return (; GPU, res, Ri, s, ν_h, ν_v, advection_scheme, horizontal_hyperviscosity, short_duration)
end
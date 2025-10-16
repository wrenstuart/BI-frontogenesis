function sim_params() 
    Ri = 1
    s = 1e4
    ν_v = 8e-3
    ν_h = 3.2e+2
    GPU = true
    res = (512, 512, 64)
    advection_scheme = CenteredSecondOrder
    horizontal_hyperviscosity = false
    short_duration = false
    return (; GPU, res, Ri, s, ν_h, ν_v, advection_scheme, horizontal_hyperviscosity, short_duration)
end
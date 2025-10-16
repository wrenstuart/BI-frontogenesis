function sim_params() 
    Ri = 1
    s = 1e4
    ν_v = 2.5e-4
    ν_h = 1e+1
    GPU = true
    res = (512, 512, 64)
    advection_scheme = CenteredSecondOrder
    horizontal_hyperviscosity = false
    short_duration = false
    return (; GPU, res, Ri, s, ν_h, ν_v, advection_scheme, horizontal_hyperviscosity, short_duration)
end
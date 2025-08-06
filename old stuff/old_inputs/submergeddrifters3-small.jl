function sim_params()
    Ri = 1
    s = 1e4
    ν_v = 1e-3
    ν_h = 1e+1
    GPU = false
    res = (64, 64, 8)
    advection_scheme = CenteredSecondOrder
    horizontal_hyperviscosity = false
    short_duration = true
    fix_drifters_below_surface = true
    return (; GPU, res, Ri, s, ν_h, ν_v, advection_scheme, horizontal_hyperviscosity, short_duration, fix_drifters_below_surface)
end
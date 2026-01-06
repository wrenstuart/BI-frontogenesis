# Now outputting halo points, which we would like to check
# This was before updating to Oceananigans v0.95
function sim_params() 
    Ri = 1
    s = 1e4
    ν_v = 1e-3
    ν_h = 4e+1
    GPU = true
    res = (768, 768, 64)
    advection_scheme = CenteredSecondOrder
    horizontal_hyperviscosity = false
    short_duration = true
    diffusive_cfl = 0.3
    return (; GPU, res, Ri, s, ν_h, ν_v, advection_scheme, horizontal_hyperviscosity, short_duration, diffusive_cfl)
end
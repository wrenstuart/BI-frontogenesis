# Some global numerical instability after running---likely need to decrease CFL number
function sim_params() 
    Ri = 1
    s = 1e4
    ν_v = 1e-3
    ν_h = 4e+1
    GPU = true
    res = (256, 256, 32)
    advection_scheme = Centered
    horizontal_hyperviscosity = false
    short_duration = false
    diffusive_cfl = 0.2
    return (; GPU, res, Ri, s, ν_h, ν_v, advection_scheme, horizontal_hyperviscosity, short_duration, diffusive_cfl)
end
# For generating strings to label simulations based on their parameters

function get_label(log_Ri, log_s)
    return "_Ri_" * string(log_Ri) * "_s_" * string(log_s)
end

function get_scales(Ri, s)

    # Ri = N²f²/M⁴
    # s = N²/f²

    # Set some dimensional parameters
    H = 50      # depth of the mixed layer (WLOG)
    f = 1e-4    # Coriolis parameter (WLOG)

    # Calculate the other non-dimensional parameters
    α = (s * Ri) ^ 0.5  # N²/M²
    s = α^2 / Ri        # N²/f²
    λ = s / α           # M²/f²
    Ro₀ = Ri ^ -0.5     # The initial Rossby number

    # Calculate the other dimensional parameters and scales
    M² = λ * f^2                # Horizontal buoyancy gradient
    N² = s * f^2                # Vertical buoyancy gradient
    L = λ * H * (1+Ri) ^ 0.5    # Horizontal lengthscale (Rossby deformation radius)
    U = λ * f * H               # Velocity scale
    T = (1+Ri) ^ 0.5 / f        # Timescale of growth (and expected timescale of later flow, if Ro ∼ Ri^(-0.5)

    dim_params = (M² = M², N² = N², L = L, U = U, T = T, f = f, H = H)
    return dim_params

end
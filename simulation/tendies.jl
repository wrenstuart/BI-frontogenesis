using Oceananigans.Biogeochemistry: biogeochemical_transition, biogeochemical_drift_velocity
using Oceananigans.TurbulenceClosures: ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ, ∂ⱼ_τ₃ⱼ, ∇_dot_qᶜ
using Oceananigans.TurbulenceClosures: immersed_∂ⱼ_τ₁ⱼ, immersed_∂ⱼ_τ₂ⱼ, immersed_∂ⱼ_τ₃ⱼ, immersed_∇_dot_qᶜ
using Oceananigans.Forcings: with_advective_forcing
using Oceananigans.Advection
using Oceananigans.Coriolis
using Oceananigans.Utils: SumOfArrays

"return the ``x and y``-gradient of hydrostatic pressure"
hydrostatic_pressure_gradient_x(i, j, k, grid, hydrostatic_pressure) = ∂xᶠᶜᶜ(i, j, k, grid, hydrostatic_pressure)
hydrostatic_pressure_gradient_y(i, j, k, grid, hydrostatic_pressure) = ∂yᶜᶠᶜ(i, j, k, grid, hydrostatic_pressure)
hydrostatic_pressure_gradient_x(i, j, k, grid, ::Nothing) = zero(grid)
hydrostatic_pressure_gradient_y(i, j, k, grid, ::Nothing) = zero(grid)

"return the ``x and y``-gradient of nonhydrostatic pressure"
nonhydrostatic_pressure_gradient_x(i, j, k, grid, nonhydrostatic_pressure) = ∂xᶠᶜᶜ(i, j, k, grid, nonhydrostatic_pressure)
nonhydrostatic_pressure_gradient_y(i, j, k, grid, nonhydrostatic_pressure) = ∂yᶜᶠᶜ(i, j, k, grid, nonhydrostatic_pressure)
nonhydrostatic_pressure_gradient_x(i, j, k, grid, ::Nothing) = zero(grid)
nonhydrostatic_pressure_gradient_y(i, j, k, grid, ::Nothing) = zero(grid)

array_to_function(arr) = (i, j, k) -> arr[i, j, k]
a2f(arr) = array_to_function(arr)

# Taking a function taking indices as inputs (rather than an array)
@inline f̅ᶜᵃᵃ(f::Function) = (i, j, k) -> (f(i+1, j, k) + f(i, j, k))/2
@inline f̅ᶠᵃᵃ(f::Function) = (i, j, k) -> (f(i-1, j, k) + f(i, j, k))/2
@inline f̅ᵃᶜᵃ(f::Function) = (i, j, k) -> (f(i, j+1, k) + f(i, j, k))/2
@inline f̅ᵃᶠᵃ(f::Function) = (i, j, k) -> (f(i, j-1, k) + f(i, j, k))/2
@inline f̅ᵃᵃᶜ(f::Function) = (i, j, k) -> (f(i, j, k+1) + f(i, j, k))/2
@inline f̅ᵃᵃᶠ(f::Function) = (i, j, k) -> (f(i, j, k-1) + f(i, j, k))/2
@inline f̅ᶜᵃᵃ(arr) = f̅ᶜᵃᵃ(a2f(arr))
@inline f̅ᶠᵃᵃ(arr) = f̅ᶠᵃᵃ(a2f(arr))
@inline f̅ᵃᶜᵃ(arr) = f̅ᵃᶜᵃ(a2f(arr))
@inline f̅ᵃᶠᵃ(arr) = f̅ᵃᶠᵃ(a2f(arr))
@inline f̅ᵃᵃᶜ(arr) = f̅ᵃᵃᶜ(a2f(arr))
@inline f̅ᵃᵃᶠ(arr) = f̅ᵃᵃᶠ(a2f(arr))

@inline δxᶜᵃᵃ(f::Function) = (i, j, k) ->  f(i+1, j, k) - f(i, j, k)
@inline δxᶠᵃᵃ(f::Function) = (i, j, k) -> -f(i-1, j, k) + f(i, j, k)
@inline δyᵃᶜᵃ(f::Function) = (i, j, k) ->  f(i, j+1, k) - f(i, j, k)
@inline δyᵃᶠᵃ(f::Function) = (i, j, k) -> -f(i, j-1, k) + f(i, j, k)
@inline δzᵃᵃᶜ(f::Function) = (i, j, k) ->  f(i, j, k+1) - f(i, j, k)
@inline δzᵃᵃᶠ(f::Function) = (i, j, k) -> -f(i, j, k-1) + f(i, j, k)
@inline δxᶜᵃᵃ(arr::Array) = δxᶜᵃᵃ(a2f(arr))
@inline δxᶠᵃᵃ(arr::Array) = δxᶠᵃᵃ(a2f(arr))
@inline δyᵃᶜᵃ(arr::Array) = δyᵃᶜᵃ(a2f(arr))
@inline δyᵃᶠᵃ(arr::Array) = δyᵃᶠᵃ(a2f(arr))
@inline δzᵃᵃᶜ(arr::Array) = δzᵃᵃᶜ(a2f(arr))
@inline δzᵃᵃᶠ(arr::Array) = δzᵃᵃᶠ(a2f(arr))

@inline f̅ᶜᶜᵃ(f) = f̅ᶜᵃᵃ(f̅ᵃᶜᵃ(f))
@inline f̅ᶜᶠᵃ(f) = f̅ᶜᵃᵃ(f̅ᵃᶠᵃ(f))
@inline f̅ᶠᶜᵃ(f) = f̅ᶠᵃᵃ(f̅ᵃᶜᵃ(f))
@inline f̅ᶠᶠᵃ(f) = f̅ᶠᵃᵃ(f̅ᵃᶠᵃ(f))
@inline f̅ᵃᶜᶜ(f) = f̅ᵃᶜᵃ(f̅ᵃᵃᶜ(f))
@inline f̅ᵃᶜᶠ(f) = f̅ᵃᶜᵃ(f̅ᵃᵃᶠ(f))
@inline f̅ᵃᶠᶜ(f) = f̅ᵃᶠᵃ(f̅ᵃᵃᶜ(f))
@inline f̅ᵃᶠᶠ(f) = f̅ᵃᶠᵃ(f̅ᵃᵃᶠ(f))
@inline f̅ᶜᵃᶜ(f) = f̅ᶜᵃᵃ(f̅ᵃᵃᶜ(f))
@inline f̅ᶜᵃᶠ(f) = f̅ᶜᵃᵃ(f̅ᵃᵃᶠ(f))
@inline f̅ᶠᵃᶜ(f) = f̅ᶠᵃᵃ(f̅ᵃᵃᶜ(f))
@inline f̅ᶠᵃᶠ(f) = f̅ᶠᵃᵃ(f̅ᵃᵃᶠ(f))

@inline f̅★ᵃᵃ(f) = f̅ᶜᵃᵃ(f̅ᶠᵃᵃ(f))
@inline f̅ᵃ★ᵃ(f) = f̅ᵃᶜᵃ(f̅ᵃᶠᵃ(f))
@inline f̅ᵃᵃ★(f) = f̅ᵃᵃᶜ(f̅ᵃᵃᶠ(f))

@inline ∂xᶜᶜᶜ_f(grid, arr) = (i, j, k) -> ∂xᶜᶜᶜ(i, j, k, grid, arr)
@inline ∂xᶜᶜᶠ_f(grid, arr) = (i, j, k) -> ∂xᶜᶜᶠ(i, j, k, grid, arr)
@inline ∂xᶜᶠᶜ_f(grid, arr) = (i, j, k) -> ∂xᶜᶠᶜ(i, j, k, grid, arr)
@inline ∂xᶜᶠᶠ_f(grid, arr) = (i, j, k) -> ∂xᶜᶠᶠ(i, j, k, grid, arr)
@inline ∂xᶠᶜᶜ_f(grid, arr) = (i, j, k) -> ∂xᶠᶜᶜ(i, j, k, grid, arr)
@inline ∂xᶠᶜᶠ_f(grid, arr) = (i, j, k) -> ∂xᶠᶜᶠ(i, j, k, grid, arr)
@inline ∂xᶠᶠᶜ_f(grid, arr) = (i, j, k) -> ∂xᶠᶠᶜ(i, j, k, grid, arr)
@inline ∂xᶠᶠᶠ_f(grid, arr) = (i, j, k) -> ∂xᶠᶠᶠ(i, j, k, grid, arr)
@inline ∂yᶜᶜᶜ_f(grid, arr) = (i, j, k) -> ∂yᶜᶜᶜ(i, j, k, grid, arr)
@inline ∂yᶜᶜᶠ_f(grid, arr) = (i, j, k) -> ∂yᶜᶜᶠ(i, j, k, grid, arr)
@inline ∂yᶜᶠᶜ_f(grid, arr) = (i, j, k) -> ∂yᶜᶠᶜ(i, j, k, grid, arr)
@inline ∂yᶜᶠᶠ_f(grid, arr) = (i, j, k) -> ∂yᶜᶠᶠ(i, j, k, grid, arr)
@inline ∂yᶠᶜᶜ_f(grid, arr) = (i, j, k) -> ∂yᶠᶜᶜ(i, j, k, grid, arr)
@inline ∂yᶠᶜᶠ_f(grid, arr) = (i, j, k) -> ∂yᶠᶜᶠ(i, j, k, grid, arr)
@inline ∂yᶠᶠᶜ_f(grid, arr) = (i, j, k) -> ∂yᶠᶠᶜ(i, j, k, grid, arr)
@inline ∂yᶠᶠᶠ_f(grid, arr) = (i, j, k) -> ∂yᶠᶠᶠ(i, j, k, grid, arr)
@inline ∂zᶜᶜᶜ_f(grid, arr) = (i, j, k) -> ∂zᶜᶜᶜ(i, j, k, grid, arr)
@inline ∂zᶜᶜᶠ_f(grid, arr) = (i, j, k) -> ∂zᶜᶜᶠ(i, j, k, grid, arr)
@inline ∂zᶜᶠᶜ_f(grid, arr) = (i, j, k) -> ∂zᶜᶠᶜ(i, j, k, grid, arr)
@inline ∂zᶜᶠᶠ_f(grid, arr) = (i, j, k) -> ∂zᶜᶠᶠ(i, j, k, grid, arr)
@inline ∂zᶠᶜᶜ_f(grid, arr) = (i, j, k) -> ∂zᶠᶜᶜ(i, j, k, grid, arr)
@inline ∂zᶠᶜᶠ_f(grid, arr) = (i, j, k) -> ∂zᶠᶜᶠ(i, j, k, grid, arr)
@inline ∂zᶠᶠᶜ_f(grid, arr) = (i, j, k) -> ∂zᶠᶠᶜ(i, j, k, grid, arr)
@inline ∂zᶠᶠᶠ_f(grid, arr) = (i, j, k) -> ∂zᶠᶠᶠ(i, j, k, grid, arr)

@inline ∂xᶜᵃᵃ_f(grid, f::Function) = (i, j, k) -> δxᶜᵃᵃ(f)(i, j, k) / grid.Δxᶜᵃᵃ
@inline ∂xᶠᵃᵃ_f(grid, f::Function) = (i, j, k) -> δxᶠᵃᵃ(f)(i, j, k) / grid.Δxᶠᵃᵃ
@inline ∂yᵃᶜᵃ_f(grid, f::Function) = (i, j, k) -> δyᵃᶜᵃ(f)(i, j, k) / grid.Δyᵃᶜᵃ
@inline ∂yᵃᶠᵃ_f(grid, f::Function) = (i, j, k) -> δyᵃᶠᵃ(f)(i, j, k) / grid.Δyᵃᶠᵃ
@inline ∂zᵃᵃᶜ_f(grid, f::Function) = (i, j, k) -> δzᵃᵃᶜ(f)(i, j, k) / grid.z.Δᵃᵃᶜ
@inline ∂zᵃᵃᶠ_f(grid, f::Function) = (i, j, k) -> δzᵃᵃᶠ(f)(i, j, k) / grid.z.Δᵃᵃᶠ

@inline add(f::Function, g::Function) = (i, j, k) -> f(i, j, k) + g(i, j, k)
@inline mult(f::Function, g::Function) = (i, j, k) -> f(i, j, k) * g(i, j, k)

@inline ∂x²_f(grid, f::Function) = (i, j, k) -> δxᶠᵃᵃ(δxᶜᵃᵃ(f))(i, j, k) / (grid.Δxᶜᵃᵃ * grid.Δxᶠᵃᵃ)
@inline ∂y²_f(grid, f::Function) = (i, j, k) -> δyᵃᶠᵃ(δyᵃᶜᵃ(f))(i, j, k) / (grid.Δyᵃᶜᵃ * grid.Δyᵃᶠᵃ)
@inline ∂z²_f(grid, f::Function) = (i, j, k) -> δzᵃᵃᶠ(δzᵃᵃᶜ(f))(i, j, k) / (grid.z.Δᵃᵃᶜ * grid.z.Δᵃᵃᶠ)
@inline ∇ₕ²_f(grid, f::Function) = add(∂x²_f(grid, f::Function), ∂y²_f(grid, f::Function))


@inline function u_t_func_full( # ∂u/∂t, different from the tendency G_u (which
                                # does not include nonhydrostatic pressure)
    i, j, k,
    grid,
    advection_scheme,
    coriolis,
    closure,
    buoyancy,
    background_fields,
    velocities,
    tracers,
    diffusivities,
    hydrostatic_pressure,
    nonhydrostatic_pressure)

    total_velocities = (u = SumOfArrays{2}(velocities.u, background_fields.velocities.u),
                        v = SumOfArrays{2}(velocities.v, background_fields.velocities.v),
                        w = SumOfArrays{2}(velocities.w, background_fields.velocities.w))
    model_fields = merge(velocities, tracers)
    return ( - div_𝐯u(i, j, k, grid, advection_scheme, total_velocities, velocities.u)
             - div_𝐯u(i, j, k, grid, advection_scheme, velocities, background_fields.velocities.u)  # Pretty sure can ignore this term
             - x_f_cross_U(i, j, k, grid, coriolis, velocities)
             - hydrostatic_pressure_gradient_x(i, j, k, grid, hydrostatic_pressure)
             - nonhydrostatic_pressure_gradient_x(i, j, k, grid, nonhydrostatic_pressure)
             - ∂ⱼ_τ₁ⱼ(i, j, k, grid, closure, diffusivities, clock, model_fields, buoyancy))

end

@inline function v_t_func_full( # Similar to u_t_func, different from tendency
    i, j, k,
    grid,
    advection_scheme,
    coriolis,
    closure,
    buoyancy,
    background_fields,
    velocities,
    tracers,
    diffusivities,
    hydrostatic_pressure,
    nonhydrostatic_pressure)

    total_velocities = (u = SumOfArrays{2}(velocities.u, background_fields.velocities.u),
                        v = SumOfArrays{2}(velocities.v, background_fields.velocities.v),
                        w = SumOfArrays{2}(velocities.w, background_fields.velocities.w))
    model_fields = merge(velocities, tracers)
    return ( - div_𝐯v(i, j, k, grid, advection_scheme, total_velocities, velocities.v)
             - div_𝐯v(i, j, k, grid, advection_scheme, velocities, background_fields.velocities.v)  # Pretty sure can ignore this term
             - y_f_cross_U(i, j, k, grid, coriolis, velocities)
             - hydrostatic_pressure_gradient_y(i, j, k, grid, hydrostatic_pressure)
             - nonhydrostatic_pressure_gradient_y(i, j, k, grid, nonhydrostatic_pressure)
             - ∂ⱼ_τ₂ⱼ(i, j, k, grid, closure, diffusivities, clock, model_fields, buoyancy))

end

@inline function u_t_func(i, j, k, grid, other_args)
    # For some reason the GPU can't handle taking the other args from the named
    # tuple into u_t_full_func in one step---writing it out in long form like
    # this avoids compiling problems
    a = other_args
    advection_scheme = a.advection_scheme
    coriolis = a.coriolis
    closure = a.closure
    buoyancy = a.buoyancy
    background_fields = a.background_fields
    velocities = a.velocities
    tracers = a.tracers
    diffusivities = a.diffusivities
    hydrostatic_pressure = a.hydrostatic_pressure
    nonhydrostatic_pressure = a.nonhydrostatic_pressure
    return u_t_func_full(i, j, k, grid, advection_scheme, coriolis,
              closure, buoyancy, background_fields, velocities, tracers,
              diffusivities, hydrostatic_pressure, nonhydrostatic_pressure)
end
@inline function v_t_func(i, j, k, grid, other_args)
    a = other_args
    advection_scheme = a.advection_scheme
    coriolis = a.coriolis
    closure = a.closure
    buoyancy = a.buoyancy
    background_fields = a.background_fields
    velocities = a.velocities
    tracers = a.tracers
    diffusivities = a.diffusivities
    hydrostatic_pressure = a.hydrostatic_pressure
    nonhydrostatic_pressure = a.nonhydrostatic_pressure
    return v_t_func_full(i, j, k, grid, advection_scheme, coriolis,
              closure, buoyancy, background_fields, velocities, tracers,
              diffusivities, hydrostatic_pressure, nonhydrostatic_pressure)
end

#=@inline function u_tendency_func_full(
    i, j, k,
    grid,
    advection_scheme,
    coriolis,
    closure,
    buoyancy,
    background_fields,
    velocities,
    tracers,
    diffusivities,
    hydrostatic_pressure)

    total_velocities = (u = SumOfArrays{2}(velocities.u, background_fields.velocities.u),
                        v = SumOfArrays{2}(velocities.v, background_fields.velocities.v),
                        w = SumOfArrays{2}(velocities.w, background_fields.velocities.w))
    model_fields = merge(velocities, tracers)
    return ( - div_𝐯u(i, j, k, grid, advection_scheme, total_velocities, velocities.u)
             - div_𝐯u(i, j, k, grid, advection_scheme, velocities, background_fields.velocities.u)  # Pretty sure can ignore this term
             - x_f_cross_U(i, j, k, grid, coriolis, velocities)
             - hydrostatic_pressure_gradient_x(i, j, k, grid, hydrostatic_pressure)
             - ∂ⱼ_τ₁ⱼ(i, j, k, grid, closure, diffusivities, clock, model_fields, buoyancy))

end

@inline function v_tendency_func_full(
    i, j, k,
    grid,
    advection_scheme,
    coriolis,
    closure,
    buoyancy,
    background_fields,
    velocities,
    tracers,
    diffusivities,
    hydrostatic_pressure)

    total_velocities = (u = SumOfArrays{2}(velocities.u, background_fields.velocities.u),
                        v = SumOfArrays{2}(velocities.v, background_fields.velocities.v),
                        w = SumOfArrays{2}(velocities.w, background_fields.velocities.w))
    model_fields = merge(velocities, tracers)
    return ( - div_𝐯v(i, j, k, grid, advection_scheme, total_velocities, velocities.v)
             - div_𝐯v(i, j, k, grid, advection_scheme, velocities, background_fields.velocities.v)  # Pretty sure can ignore this term
             - y_f_cross_U(i, j, k, grid, coriolis, velocities)
             - hydrostatic_pressure_gradient_y(i, j, k, grid, hydrostatic_pressure)
             - ∂ⱼ_τ₂ⱼ(i, j, k, grid, closure, diffusivities, clock, model_fields, buoyancy))

end

@inline function u_tendency_func(i, j, k, grid, other_args)
    a = other_args
    advection_scheme = a.advection_scheme
    coriolis = a.coriolis
    closure = a.closure
    buoyancy = a.buoyancy
    background_fields = a.background_fields
    velocities = a.velocities
    tracers = a.tracers
    diffusivities = a.diffusivities
    hydrostatic_pressure = a.hydrostatic_pressure
    return u_tendency_func_full(i, j, k, grid, advection_scheme, coriolis, closure, buoyancy,
                            background_fields, velocities, tracers, diffusivities, hydrostatic_pressure)
end
@inline function v_tendency_func(i, j, k, grid, other_args)
    a = other_args
    advection_scheme = a.advection_scheme
    coriolis = a.coriolis
    closure = a.closure
    buoyancy = a.buoyancy
    background_fields = a.background_fields
    velocities = a.velocities
    tracers = a.tracers
    diffusivities = a.diffusivities
    hydrostatic_pressure = a.hydrostatic_pressure
    return v_tendency_func_full(i, j, k, grid, advection_scheme, coriolis, closure, buoyancy,
                            background_fields, velocities, tracers, diffusivities, hydrostatic_pressure)
end=#

@inline function u_adv_func(i, j, k, grid, other_args)
    a = other_args
    advection_scheme = a.advection_scheme
    background_fields = a.background_fields
    velocities = a.velocities
    total_velocities = (u = SumOfArrays{2}(velocities.u, background_fields.velocities.u),
                        v = SumOfArrays{2}(velocities.v, background_fields.velocities.v),
                        w = SumOfArrays{2}(velocities.w, background_fields.velocities.w))
    return (  div_𝐯u(i, j, k, grid, advection_scheme, total_velocities, velocities.u)
            + div_𝐯u(i, j, k, grid, advection_scheme, velocities, background_fields.velocities.u))
end
@inline function u_cor_func(i, j, k, grid, other_args)
    a = other_args
    return - x_f_cross_U(i, j, k, grid, a.coriolis, a.velocities)
end
@inline function v_cor_func(i, j, k, grid, other_args)
    a = other_args
    return - y_f_cross_U(i, j, k, grid, a.coriolis, a.velocities)
end
@inline function u_visc_func(i, j, k, grid, other_args)
    a = other_args
    model_fields = merge(a.velocities, a.tracers)
    return - ∂ⱼ_τ₁ⱼ(i, j, k, grid, a.closure, a.diffusivities, clock, model_fields, a.buoyancy)
end
@inline function v_visc_func(i, j, k, grid, other_args)
    a = other_args
    model_fields = merge(a.velocities, a.tracers)
    return - ∂ⱼ_τ₂ⱼ(i, j, k, grid, a.closure, a.diffusivities, clock, model_fields, a.buoyancy)
end
@inline function u_prs_func(i, j, k, grid, other_args)
    a = other_args
    return ( - hydrostatic_pressure_gradient_x(i, j, k, grid, a.hydrostatic_pressure)
             - nonhydrostatic_pressure_gradient_x(i, j, k, grid, a.nonhydrostatic_pressure))
end
@inline function v_prs_func(i, j, k, grid, other_args)
    a = other_args
    return ( - hydrostatic_pressure_gradient_y(i, j, k, grid, a.hydrostatic_pressure)
             - nonhydrostatic_pressure_gradient_y(i, j, k, grid, a.nonhydrostatic_pressure))
end

@inline function u_err_func(i, j, k, grid, other_args)   # Error from ∇⋅(𝐮u) ≠ 𝐮⋅∇u

    a = other_args

    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w)

    return (  f̅★ᵃᵃ(u)(i, j, k) * f̅ᶠᵃᵃ(∂xᶜᶜᶜ_f(grid, u))(i, j, k)
            + f̅ᵃ★ᵃ(u)(i, j, k) * f̅ᶠᵃᵃ(∂yᶜᶜᶜ_f(grid, v))(i, j, k)
            + f̅ᵃᵃ★(u)(i, j, k) * f̅ᶠᵃᵃ(∂zᶜᶜᶜ_f(grid, w))(i, j, k))

end

@inline function v_err_func(i, j, k, grid, other_args)   # Error from ∇⋅(𝐮v) ≠ 𝐮⋅∇v

    a = other_args

    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w)

    return (  f̅★ᵃᵃ(v)(i, j, k) * f̅ᵃᶠᵃ(∂xᶜᶜᶜ_f(grid, u))(i, j, k)
            + f̅ᵃ★ᵃ(v)(i, j, k) * f̅ᵃᶠᵃ(∂yᶜᶜᶜ_f(grid, v))(i, j, k)
            + f̅ᵃᵃ★(v)(i, j, k) * f̅ᵃᶠᵃ(∂zᶜᶜᶜ_f(grid, w))(i, j, k))

end

########################
# Lagrangian ζ t-deriv #
########################

@inline function vtcl_curl_func(u_func, v_func)

    return (i, j, k, grid, other_args) -> begin
        u_f = (i, j, k) -> u_func(i, j, k, grid, other_args)
        v_f = (i, j, k) -> v_func(i, j, k, grid, other_args)
        ∂xᶠᵃᵃ_f(grid, v_f)(i, j, k) - ∂yᵃᶠᵃ_f(grid, u_f)(i, j, k)
    end

end

@inline function F_ζ_hor_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w)

    return (- f̅ᶠᶠᵃ(∂xᶜᶜᶜ_f(grid, u))(i, j, k) * f̅★ᵃᵃ(∂xᶠᶠᶜ_f(grid, v))(i, j, k)
            - f̅ᶠᶠᵃ(∂yᶜᶜᶜ_f(grid, v))(i, j, k) * f̅ᵃ★ᵃ(∂xᶠᶠᶜ_f(grid, v))(i, j, k)
            + f̅ᶠᶠᵃ(∂xᶜᶜᶜ_f(grid, u))(i, j, k) * f̅★ᵃᵃ(∂yᶠᶠᶜ_f(grid, u))(i, j, k)
            + f̅ᶠᶠᵃ(∂yᶜᶜᶜ_f(grid, v))(i, j, k) * f̅ᵃ★ᵃ(∂yᶠᶠᶜ_f(grid, u))(i, j, k))

end

@inline function F_ζ_vrt_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w)

    return (- f̅ᵃᶠᶜ(∂xᶠᶜᶠ_f(grid, w))(i, j, k) * f̅ᶠᵃᶜ(∂zᶜᶠᶠ_f(grid, v))(i, j, k)
            + f̅ᶠᵃᶜ(∂yᶜᶠᶠ_f(grid, w))(i, j, k) * f̅ᵃᶠᶜ(∂zᶠᶜᶠ_f(grid, u))(i, j, k))
    
end

@inline function ζ_adv_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w)
    ζ_f = (i, j, k) -> ∂xᶠᶠᶜ(i, j, k, grid, v) - ∂yᶠᶠᶜ(i, j, k, grid, u)

    return (  f̅★ᵃᵃ(f̅ᵃᶠᵃ(u))(i, j, k) * f̅ᶠᵃᵃ(∂xᶜᵃᵃ_f(grid, ζ_f))(i, j, k)
            + f̅ᵃ★ᵃ(f̅ᶠᵃᵃ(v))(i, j, k) * f̅ᵃᶠᵃ(∂yᵃᶜᵃ_f(grid, ζ_f))(i, j, k)
            + f̅ᶠᶠᵃ(f̅ᵃᵃᶜ(w))(i, j, k) * f̅ᵃᵃᶜ(∂zᵃᵃᶠ_f(grid, ζ_f))(i, j, k))
    
end

@inline function F_ζ_cor_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    f = a.coriolis.f

    δ_f = (i, j, k) -> ∂xᶜᶜᶜ(i, j, k, grid, u) + ∂yᶜᶜᶜ(i, j, k, grid, v)
    return - f * f̅ᶠᶠᵃ(δ_f)(i, j, k)
    
end

@inline function ζ_h_visc_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    ζ_f = (i, j, k) -> ∂xᶠᶠᶜ(i, j, k, grid, v) - ∂yᶠᶠᶜ(i, j, k, grid, u)

    return other_args.diffusivities[1].ν * ∇ₕ²_f(grid, ζ_f)(i, j, k)
    
end

@inline function ζ_v_visc_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    ζ_f = (i, j, k) -> ∂xᶠᶠᶜ(i, j, k, grid, v) - ∂yᶠᶠᶜ(i, j, k, grid, u)

    return other_args.diffusivities[2].ν * ∂z²_f(grid, ζ_f)(i, j, k)
    
end

@inline ζ_err_func = vtcl_curl_func(u_err_func, v_err_func)

@inline ζ_t_func = vtcl_curl_func(u_t_func, v_t_func)

########################
# Lagrangian δ t-deriv #
########################

@inline function hor_div_func(u_func, v_func)

    return (i, j, k, grid, other_args) -> begin
        u_f = (i, j, k) -> u_func(i, j, k, grid, other_args)
        v_f = (i, j, k) -> v_func(i, j, k, grid, other_args)
        ∂xᶜᵃᵃ_f(grid, u_f)(i, j, k) + ∂yᵃᶜᵃ_f(grid, v_f)(i, j, k)
    end

end

@inline function u_func(i, j, k)
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    return u[i, j, k]
end
@inline function v_func(i, j, k)
    v = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.v)
    return v[i, j, k]
end
@inline δ_func = hor_div_func(u_func, v_func)

@inline δ_t_func     = hor_div_func(u_t_func, v_t_func)
@inline δ_err_func   = hor_div_func(u_err_func, v_err_func)
@inline F_δ_cor_func = hor_div_func(u_cor_func, v_cor_func)
@inline F_δ_prs_func = hor_div_func(u_prs_func, v_prs_func)

@inline function δ_h_visc_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    δ_f = (i, j, k) -> ∂xᶜᶜᶜ(i, j, k, grid, u) + ∂yᶜᶜᶜ(i, j, k, grid, v)

    return other_args.diffusivities[1].ν * ∇ₕ²_f(grid, δ_f)(i, j, k)
    
end

@inline function δ_v_visc_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    δ_f = (i, j, k) -> ∂xᶜᶜᶜ(i, j, k, grid, u) + ∂yᶜᶜᶜ(i, j, k, grid, v)

    return other_args.diffusivities[2].ν * ∂z²_f(grid, δ_f)(i, j, k)
    
end

@inline function δ_adv_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w)
    u_xx_f = ∂xᶠᵃᵃ_f(grid, ∂xᶜᶜᶜ_f(grid, u))
    u_yx_f = ∂yᵃᶠᵃ_f(grid, ∂xᶜᶜᶜ_f(grid, u))
    u_zx_f = ∂zᵃᵃᶠ_f(grid, ∂xᶜᶜᶜ_f(grid, u))
    v_xy_f = ∂xᶠᵃᵃ_f(grid, ∂yᶜᶜᶜ_f(grid, v))
    v_yy_f = ∂yᵃᶠᵃ_f(grid, ∂yᶜᶜᶜ_f(grid, v))
    v_zy_f = ∂zᵃᵃᶠ_f(grid, ∂yᶜᶜᶜ_f(grid, v))

    return (  f̅★ᵃᵃ(f̅ᶜᵃᵃ(u))(i, j, k) * f̅ᶜᵃᵃ(u_xx_f)(i, j, k)
            + f̅★ᵃᵃ(f̅ᵃᶜᵃ(v))(i, j, k) * f̅ᵃᶜᵃ(u_yx_f)(i, j, k)
            + f̅★ᵃᵃ(f̅ᵃᵃᶜ(w))(i, j, k) * f̅ᵃᵃᶜ(u_zx_f)(i, j, k)
            + f̅ᵃ★ᵃ(f̅ᶜᵃᵃ(u))(i, j, k) * f̅ᶜᵃᵃ(v_xy_f)(i, j, k)
            + f̅ᵃ★ᵃ(f̅ᵃᶜᵃ(v))(i, j, k) * f̅ᵃᶜᵃ(v_yy_f)(i, j, k)
            + f̅ᵃ★ᵃ(f̅ᵃᵃᶜ(w))(i, j, k) * f̅ᵃᵃᶜ(v_zy_f)(i, j, k))

end

@inline function F_δ_hor_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)

    return - (  (∂xᶜᵃᵃ_f(grid, f̅★ᵃᵃ(u)))(i, j, k) ^ 2
              + 2∂xᶜᵃᵃ_f(grid, f̅ᶠᶜᵃ(v))(i, j, k) * ∂yᵃᶜᵃ_f(grid, f̅ᶜᶠᵃ(u))(i, j, k)
              + (∂yᵃᶜᵃ_f(grid, f̅ᵃ★ᵃ(v)))(i, j, k) ^ 2)
    
end

@inline function F_δ_vrt_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w)

    return - (  ∂xᶜᵃᵃ_f(grid, f̅ᶠᵃᶜ(w))(i, j, k) * ∂zᵃᵃᶜ_f(grid, f̅ᶜᵃᶠ(u))(i, j, k)
              + ∂yᵃᶜᵃ_f(grid, f̅ᵃᶠᶜ(w))(i, j, k) * ∂zᵃᵃᶜ_f(grid, f̅ᵃᶜᶠ(v))(i, j, k))

end

#=@inline function u_div𝐯_func(i, j, k, grid, other_args)
    
    a = other_args
    total_velocities = (u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u),
                        v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v),
                        w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w))
    return (  div_𝐯u(i, j, k, grid, a.advection_scheme, total_velocities, a.velocities.u)
            + div_𝐯u(i, j, k, grid, a.advection_scheme, a.velocities, a.background_fields.velocities.u))
    
end

@inline function v_div𝐯_func(i, j, k, grid, other_args)
    
    a = other_args
    total_velocities = (u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u),
                        v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v),
                        w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w))
    return (  div_𝐯v(i, j, k, grid, a.advection_scheme, total_velocities, a.velocities.v)
            + div_𝐯v(i, j, k, grid, a.advection_scheme, a.velocities, a.background_fields.velocities.v))
    
end

@inline function my_u_div𝐯_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w)
    return (  ∂xᶠᵃᵃ_f(grid, mult(f̅ᶜᵃᵃ(u), f̅ᶜᵃᵃ(u)))(i, j, k)
            + ∂yᵃᶜᵃ_f(grid, mult(f̅ᶠᵃᵃ(v), f̅ᵃᶠᵃ(u)))(i, j, k)
            + ∂zᵃᵃᶜ_f(grid, mult(f̅ᶠᵃᵃ(w), f̅ᵃᵃᶠ(u)))(i, j, k))
    
end

@inline function my_v_div𝐯_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    w = SumOfArrays{2}(a.velocities.w, a.background_fields.velocities.w)
    return (  ∂xᶜᵃᵃ_f(grid, mult(f̅ᵃᶠᵃ(u), f̅ᶠᵃᵃ(v)))(i, j, k)
            + ∂yᵃᶠᵃ_f(grid, mult(f̅ᵃᶜᵃ(v), f̅ᵃᶜᵃ(v)))(i, j, k)
            + ∂zᵃᵃᶜ_f(grid, mult(f̅ᵃᶠᵃ(w), f̅ᵃᵃᶠ(v)))(i, j, k))
    
end=#

#=@inline function ζ_h_adv_func(i, j, k, grid, other_args)

    a = other_args
    u = SumOfArrays{2}(a.velocities.u, a.background_fields.velocities.u)
    v = SumOfArrays{2}(a.velocities.v, a.background_fields.velocities.v)
    ζ_f = (i, j, k) -> ∂xᶠᶠᶜ(i, j, k, grid, v) - ∂yᶠᶠᶜ(i, j, k, grid, u)

    return (  f̅★ᵃᵃ(f̅ᵃᶠᵃ(u))(i, j, k) * f̅ᶠᵃᵃ(∂xᶜᵃᵃ_f(grid, ζ_f))(i, j, k)
            + f̅ᵃ★ᵃ(f̅ᶠᵃᵃ(v))(i, j, k) * f̅ᵃᶠᵃ(∂yᵃᶜᵃ_f(grid, ζ_f))(i, j, k))
    
end=#
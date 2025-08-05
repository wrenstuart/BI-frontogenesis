# For diagnosing velocity gradient budgets
# (in particular, the vorticity equation)

# It turns out that the budget holds perfectly
# at Eulerian gridpoints!!!

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.BoundaryConditions
using Oceananigans.Operators
using Oceananigans.Fields
using CairoMakie
using CUDA
include("simulation/tendies2.jl")

Nz = 10
grid = RectilinearGrid(CPU(), topology = (Periodic, Periodic, Bounded), size = (Nz, Nz, Nz), x = (0, 2π), y = (0, 2π), z = (0, 2π))

i = 2
j = 3
k = 5
x = grid.xᶜᵃᵃ[i]
y = grid.yᵃᶜᵃ[j]
z = grid.zᵃᵃᶜ[k]
# Pick a point

# Next try it with full fields? KernelFunctionOperation or whatever

velocities = VelocityFields(grid)
u, v, w = velocities

#=

# Version where I calculated ∇⋅(𝐮u) both analytically and using Oceananigans discretisation

k = 1
α = 0.1
u_back_func(x, y, z, t) = α*z
set!(u, (x, y, z) -> cos(k*x)cos(k*z))
set!(v, (x, y, z) -> 0)
set!(w, (x, y, z) -> -sin(k*x)sin(k*z))
fill_halo_regions!(u, v, w)
u_div𝐯_func(x, y, z) = -k * (sin(2k*x)cos(k*z)^2 + 2sin(k*x)cos(k*z)*α*z + sin(k*x)cos(k*x)cos(2k*z) + α*z*sin(k*x)cos(k*z)) - α*sin(k*x)sin(k*z)=#

α = 0.1
u_back_func(x, y, z, t) = α*z
set!(u, (x, y, z) -> sin(x)cos(y)cos(z) + cos(x)sin(y)cos(z))
set!(v, (x, y, z) -> -cos(x)sin(y)cos(z))
set!(w, (x, y, z) -> sin(x)sin(y)sin(z))
fill_halo_regions!(u, v, w)
u_div𝐯_func(x, y, z) = 0

tracers = TracerFields((:b,), grid)
pHY′ = CenterField(grid)
# pNHS = CenterField(grid)
advection_scheme = CenteredSecondOrder()
#diff_h = ScalarBiharmonicDiffusivity(Oceananigans.TurbulenceClosures.HorizontalFormulation(), Float64, ν = 1e1, κ = 1e1)
#diff_v = ScalarBiharmonicDiffusivity(Oceananigans.TurbulenceClosures.VerticalFormulation(), Float64, ν = 1e-3, κ = 1e-3)
diff_h = ScalarBiharmonicDiffusivity(Oceananigans.TurbulenceClosures.HorizontalFormulation(), Float64, ν = 0.1, κ = 0.1)
diff_v = ScalarBiharmonicDiffusivity(Oceananigans.TurbulenceClosures.VerticalFormulation(), Float64, ν = 0.1, κ = 0.1)
closure = (diff_h, diff_v)
#diffusivities = ((ν = 1e1, κ = 1e1), (ν = 1e-3, κ = 1e-3))
diffusivities = ((ν = 0.1, κ = 0.1), (ν = 0.1, κ = 0.1))
u_back = Field{Face, Center, Center}(grid)
v_back = Field{Center, Face, Center}(grid)
w_back = Field{Center, Center, Face}(grid)
b_back = Field{Center, Center, Center}(grid)
set!(u_back, (x, y, z) -> u_back_func(x, y, z, 0))
fill_halo_regions!(u_back)
background_fields = (velocities = (u = u_back, v = v_back, w = w_back),
                     tracers = (b_back))
other_args = (advection_scheme = advection_scheme,
              coriolis = FPlane(f = 1e+0),
              closure = closure,
              buoyancy = Buoyancy(model = BuoyancyTracer()),
              background_fields = background_fields,
              velocities = velocities,
              tracers = tracers,
              diffusivities = diffusivities,
              hydrostatic_pressure = pHY′)

@inline u_tendency_op = KernelFunctionOperation{Face, Center, Center}(u_tendency_func, grid, other_args)
@inline v_tendency_op = KernelFunctionOperation{Center, Face, Center}(v_tendency_func, grid, other_args)
@inline u_visc_op     = KernelFunctionOperation{Face, Center, Center}(u_visc_func,     grid, other_args)
@inline v_visc_op     = KernelFunctionOperation{Center, Face, Center}(v_visc_func,     grid, other_args)
@inline u_cor_op      = KernelFunctionOperation{Face, Center, Center}(u_cor_func,      grid, other_args)
@inline v_cor_op      = KernelFunctionOperation{Center, Face, Center}(v_cor_func,      grid, other_args)
@inline u_div𝐯_op     = KernelFunctionOperation{Face, Center, Center}(   u_div𝐯_func, grid, other_args)
@inline v_div𝐯_op     = KernelFunctionOperation{Center, Face, Center}(   v_div𝐯_func, grid, other_args)
@inline my_u_div𝐯_op  = KernelFunctionOperation{Face, Center, Center}(my_u_div𝐯_func, grid, other_args)
@inline my_v_div𝐯_op  = KernelFunctionOperation{Center, Face, Center}(my_v_div𝐯_func, grid, other_args)
@inline ζ_adv_op      = KernelFunctionOperation{Face,   Face, Center}(    ζ_adv_func, grid, other_args)
@inline F_ζ_hor_op    = KernelFunctionOperation{Face,   Face, Center}(  F_ζ_hor_func, grid, other_args)
@inline F_ζ_vrt_op    = KernelFunctionOperation{Face,   Face, Center}(  F_ζ_vrt_func, grid, other_args)
@inline u_err_op      = KernelFunctionOperation{Face, Center, Center}(    u_err_func, grid, other_args)
@inline v_err_op      = KernelFunctionOperation{Center, Face, Center}(    v_err_func, grid, other_args)
u_tendency = Field(u_tendency_op)
v_tendency = Field(v_tendency_op)
u_visc     = Field(u_visc_op)
v_visc     = Field(v_visc_op)
u_cor      = Field(u_cor_op)
v_cor      = Field(v_cor_op)
u_div𝐯     = Field(u_div𝐯_op)
v_div𝐯     = Field(v_div𝐯_op)
my_u_div𝐯  = Field(my_u_div𝐯_op)
my_v_div𝐯  = Field(my_v_div𝐯_op)
ζ_adv      = Field(ζ_adv_op)
F_ζ_hor    = Field(F_ζ_hor_op)
F_ζ_vrt    = Field(F_ζ_vrt_op)
u_err      = Field(u_err_op)
v_err      = Field(v_err_op)
compute!(u_tendency)
compute!(v_tendency)
compute!(u_visc)
compute!(v_visc)
compute!(u_cor)
compute!(v_cor)
compute!(u_div𝐯)
compute!(v_div𝐯)
compute!(my_u_div𝐯)
compute!(my_v_div𝐯)
compute!(ζ_adv)
compute!(F_ζ_hor)
compute!(F_ζ_vrt)
compute!(u_err)
compute!(v_err)
ζ = ∂x(v) - ∂y(u)
ζ_div𝐯 = ∂x(v_div𝐯) - ∂y(u_div𝐯)
ζ_err  = ∂x(v_err)  - ∂y(u_err)
ζ_tendency = ∂x(v_tendency) - ∂y(u_tendency)
ζ_visc = ∂x(v_visc) - ∂y(u_visc)
ζ_cor = ∂x(v_cor) - ∂y(u_cor)
nothing

# Get plottable arrays from Oceananigans fields

zs = parent(grid.zᵃᵃᶜ)
ζ_tendency_func(z₀) = interpolate((x, y, z₀), ζ_tendency)
ζ_visc_func(z₀) = interpolate((x, y, z₀), ζ_visc)
ζ_cor_func(z₀) = interpolate((x, y, z₀), ζ_cor)
ζ_div𝐯_func(z₀) = interpolate((x, y, z₀), ζ_div𝐯)
ζ_adv_func(z₀) = interpolate((x, y, z₀), ζ_adv)
F_ζ_hor_func(z₀) = interpolate((x, y, z₀), F_ζ_hor)
F_ζ_vrt_func(z₀) = interpolate((x, y, z₀), F_ζ_vrt)
ζ_err_func(z₀) = interpolate((x, y, z₀), ζ_err)
ζ_tendencys = ζ_tendency_func.(zs)
ζ_viscs = ζ_visc_func.(zs)
ζ_cors = ζ_cor_func.(zs)
ζ_div𝐯s = ζ_div𝐯_func.(zs)
ζ_advs = ζ_adv_func.(zs)
F_ζ_hors = F_ζ_hor_func.(zs)
F_ζ_vrts = F_ζ_vrt_func.(zs)
ζ_errs = ζ_err_func.(zs)

fig = Figure()
ax = Axis(fig[1, 1], title = "𝐳̂⋅∇×(∇⋅(𝐮𝐮))")
# lines!(ax, zs, ζ_div𝐯s, label = "𝐳̂⋅∇×(∇⋅(𝐮𝐮)) (theirs)")
# lines!(ax, zs, (ζ_advs - F_ζ_hors - F_ζ_vrts + ζ_errs), label = "𝐳̂⋅∇×(∇⋅(𝐮𝐮)) (mine)")
lines!(ax, zs, ζ_tendencys, label = "ζ_tendency (theirs)")
lines!(ax, zs, -ζ_advs - ζ_errs + F_ζ_hors + F_ζ_vrts + ζ_viscs + ζ_cors, label = "ζ_tendency (mine)")
#=lines!(ax, zs, ζ_advs, label = "ζ_adv", linestyle = :dash)
lines!(ax, zs, F_ζ_hors, label = "F_ζ_hor", linestyle = :dash)
lines!(ax, zs, F_ζ_vrts, label = "F_ζ_vrt", linestyle = :dash)
lines!(ax, zs, ζ_errs, label = "ζ_err", linestyle = :dot)
lines!(ax, zs, ζ_tendencys, label = "ζ_tendency", linestyle = :dot)
lines!(ax, zs, ζ_viscs, label = "ζ_visc", linestyle = :dot)
lines!(ax, zs, ζ_cors, label = "ζ_cor", linestyle = :dot)=#
Legend(fig[1, 2], ax, title = "Legend", position = :bottomright, fontsize = 10, colorbar = false)
display(fig)

# Note that the budget doesn't hold at i ≤ 1 in this worked example,
# but it seems to hold fine for a running simulation. Not sure why.

#=diffs = zeros(Float64, 20)
x_p, y_p, z_p = rand(Float64, 3) * 2π
ζ_tendency_p = interpolate((x_p, y_p, z_p), ζ_tendency)
ζ_adv_p = interpolate((x_p, y_p, z_p), ζ_adv)
ζ_err_p = interpolate((x_p, y_p, z_p), ζ_err)
F_ζ_hor_p = interpolate((x_p, y_p, z_p), F_ζ_hor)
F_ζ_vrt_p = interpolate((x_p, y_p, z_p), F_ζ_vrt)
ζ_cor_p = interpolate((x_p, y_p, z_p), ζ_cor)
ζ_visc_p = interpolate((x_p, y_p, z_p), ζ_visc)
fields = (; ζ_tendency, ζ_adv, ζ_err, F_ζ_hor, F_ζ_vrt, ζ_cor, ζ_visc)
@info ζ_tendency_p + ζ_adv_p + ζ_err_p - F_ζ_hor_p - F_ζ_vrt_p - ζ_cor_p - ζ_visc_p
for key in keys(fields)
    field = fields[key]
    if field isa Oceananigans.AbstractOperations.BinaryOperation
        @info key
    end
end
# Works perfectly=#

function grid_interpolateᶠᶠᶜ(grid, field, x::Float64, y::Float64)   # Interpolate var to surface position (x, y) between gridpoints


    i₋ = Int(floor(x/grid.Lx * grid.Nx)) + 1
    j₋ = Int(floor(y/grid.Ly * grid.Ny)) + 1
    @info i₋, j₋
    i₊ = i₋ % grid.Nx + 1
    j₊ = j₋ % grid.Ny + 1
    x_frac = (x - grid.xᶠᵃᵃ[i₋]) / grid.Δxᶠᵃᵃ
    y_frac = (y - grid.yᵃᶠᵃ[j₋]) / grid.Δyᵃᶠᵃ
    @info x_frac, y_frac

    f₋₋ = field[i₋, j₋, 1]
    f₋₊ = field[i₋, j₊, 1]
    f₊₋ = field[i₊, j₋, 1]
    f₊₊ = field[i₊, j₊, 1]
    f = (1-x_frac) * (1-y_frac) * f₋₋ + (1-x_frac) * y_frac * f₋₊ + x_frac * (1-y_frac) * f₊₋ + x_frac * y_frac * f₊₊
    
    return f

end
# For diagnosing velocity gradient budgets
# (in particular, the vorticity equation)

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
diff_h = ScalarBiharmonicDiffusivity(Oceananigans.TurbulenceClosures.HorizontalFormulation(), Float64, ν = 1e1, κ = 1e1)
diff_v = ScalarBiharmonicDiffusivity(Oceananigans.TurbulenceClosures.VerticalFormulation(), Float64, ν = 1e-3, κ = 1e-3)
closure = (diff_h, diff_v)
diffusivities = ((ν = 1e1, κ = 1e1), (ν = 1e-3, κ = 1e-3))
u_back = Field{Face, Center, Center}(grid)
v_back = Field{Center, Face, Center}(grid)
w_back = Field{Center, Center, Face}(grid)
b_back = Field{Center, Center, Center}(grid)
set!(u_back, (x, y, z) -> u_back_func(x, y, z, 0))
fill_halo_regions!(u_back)
background_fields = (velocities = (u = u_back, v = v_back, w = w_back),
                     tracers = (b_back))
other_args = (advection_scheme = advection_scheme,
              coriolis = FPlane(f = 1e-4),
              closure = closure,
              buoyancy = Buoyancy(model = BuoyancyTracer()),
              background_fields = background_fields,
              velocities = velocities,
              tracers = tracers,
              diffusivities = diffusivities,
              hydrostatic_pressure = pHY′)

@inline u_div𝐯_op     = KernelFunctionOperation{Face, Center, Center}(   u_div𝐯_func, grid, other_args)
@inline v_div𝐯_op     = KernelFunctionOperation{Center, Face, Center}(   v_div𝐯_func, grid, other_args)
@inline my_u_div𝐯_op  = KernelFunctionOperation{Face, Center, Center}(my_u_div𝐯_func, grid, other_args)
@inline my_v_div𝐯_op  = KernelFunctionOperation{Center, Face, Center}(my_v_div𝐯_func, grid, other_args)
@inline ζ_adv_op      = KernelFunctionOperation{Face,   Face, Center}(    ζ_adv_func, grid, other_args)
@inline F_ζ_hor_op    = KernelFunctionOperation{Face,   Face, Center}(  F_ζ_hor_func, grid, other_args)
@inline F_ζ_vrt_op    = KernelFunctionOperation{Face,   Face, Center}(  F_ζ_vrt_func, grid, other_args)
@inline u_err_op      = KernelFunctionOperation{Face, Center, Center}(    u_err_func, grid, other_args)
@inline v_err_op      = KernelFunctionOperation{Center, Face, Center}(    v_err_func, grid, other_args)
u_div𝐯     = Field(u_div𝐯_op)
v_div𝐯     = Field(v_div𝐯_op)
my_u_div𝐯  = Field(my_u_div𝐯_op)
my_v_div𝐯  = Field(my_v_div𝐯_op)
ζ_adv      = Field(ζ_adv_op)
F_ζ_hor    = Field(F_ζ_hor_op)
F_ζ_vrt    = Field(F_ζ_vrt_op)
u_err      = Field(u_err_op)
v_err      = Field(v_err_op)
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
F_ζ_vrt_2 = - ∂x(w) * ∂z(v) + ∂y(w) * ∂z(u)
F_ζ_hor_2 = - (∂x(u) + ∂y(v)) * ζ
ζ_adv_2 = u * ∂x(ζ) + v * ∂y(ζ) + w * ∂z(ζ)
nothing

# Get plottable arrays from Oceananigans fields

zs = parent(grid.zᵃᵃᶜ)
ζ_div𝐯_func(z₀) = interpolate((x, y, z₀), ζ_div𝐯)
ζ_adv_func(z₀) = interpolate((x, y, z₀), ζ_adv)
F_ζ_hor_func(z₀) = interpolate((x, y, z₀), F_ζ_hor)
F_ζ_vrt_func(z₀) = interpolate((x, y, z₀), F_ζ_vrt)
ζ_err_func(z₀) = interpolate((x, y, z₀), ζ_err)
ζ_div𝐯s = ζ_div𝐯_func.(zs)
ζ_advs = ζ_adv_func.(zs)
F_ζ_hors = F_ζ_hor_func.(zs)
F_ζ_vrts = F_ζ_vrt_func.(zs)
ζ_errs = ζ_err_func.(zs)

fig = Figure()
ax = Axis(fig[1, 1], title = "𝐳̂⋅∇×(∇⋅(𝐮𝐮))")
lines!(ax, zs, ζ_div𝐯s, label = "𝐳̂⋅∇×(∇⋅(𝐮𝐮)) (theirs)")
lines!(ax, zs, (ζ_advs - F_ζ_hors - F_ζ_vrts + ζ_errs), label = "𝐳̂⋅∇×(∇⋅(𝐮𝐮)) (mine)")
lines!(ax, zs, ζ_advs, label = "ζ_adv", linestyle = :dash)
lines!(ax, zs, F_ζ_hors, label = "F_ζ_hor", linestyle = :dash)
lines!(ax, zs, F_ζ_vrts, label = "F_ζ_vrt", linestyle = :dash)
lines!(ax, zs, ζ_errs, label = "ζ_err", linestyle = :dot)
Legend(fig[1, 2], ax, title = "Legend", position = :bottomright, fontsize = 10, colorbar = false)
display(fig)

# Note that the budget doesn't hold at i ≤ 1 in this worked example,
# but it seems to hold fine for a running simulation. Not sure why.
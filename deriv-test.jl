using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.BoundaryConditions
using Oceananigans.Operators
using CairoMakie
using CUDA

#Oceananigans.Operators: ∂z²???

Nz = 20
grid_gpu = RectilinearGrid(GPU(), topology = (Periodic, Periodic, Bounded), size = (10, 10, Nz), x = (0, 10), y = (0, 10), z = (0, 10))
grid_cpu = RectilinearGrid(CPU(), topology = (Periodic, Periodic, Bounded), size = (10, 10, Nz), x = (0, 10), y = (0, 10), z = (0, 10))
u = Field{Face, Center, Center}(grid_gpu)
v = Field{Center, Face, Center}(grid_gpu)
#set!(u, (x, y, z) -> randn())
k = 2π/10
f(z) = cos(k*z)
set!(v, (x, y, z) -> x*f(z) + y^3*f(z)/6)
set!(u, (x, y, z) -> z)
fill_halo_regions!(u)
fill_halo_regions!(v)
ζ = ∂x(v) - ∂y(u)
δ = ∂x(u) + ∂y(v)
ζ_zz = ∂z(∂z(ζ))
#ζ_zz = ∂z²(ζ)

@inline δ²zᵃᵃᶜ(i, j, k, grid, u) = @inbounds u[i, j, k+1] - 2u[i, j, k] + u[i, j, k-1]
#@inline δ²zᵃᵃᶜ(i, j, k, grid, u) = @inbounds 1.5u[i, j, k+1] - 3.5u[i, j, k] + 2.5u[i, j, k-1]-0.5u[i, j, k-2]
@inline ∂²zᵃᵃᶜ(i, j, k, grid, u) = @inbounds δ²zᵃᵃᶜ(i, j, k, grid, u) / Δzᶠᶠᶜ(i, j, k, grid)^2
ζ_zz_op = KernelFunctionOperation{Face, Face, Center}(∂²zᵃᵃᶜ, grid_gpu, ζ)
ζ_zz_new = Field(ζ_zz_op, boundary_conditions = FieldBoundaryConditions(
    east = PeriodicBoundaryCondition(),
    west = PeriodicBoundaryCondition(),
    north = PeriodicBoundaryCondition(),
    south = PeriodicBoundaryCondition(),
    top = GradientBoundaryCondition(0),
    bottom = GradientBoundaryCondition(0)))
compute!(ζ_zz_new)

@inline function ∂²zᵃᵃᶠ_top(i, j, k, grid, u)
    Δz² = Δzᶠᶠᶜ(i, j, k, grid)^2
    δ²z = (3u[i, j, k] - 7u[i, j, k-1] + 5u[i, j, k-2] - u[i, j, k-3]) / 2
    return δ²z / Δz²
end
ζ_zz_op_face = KernelFunctionOperation{Face, Face, Face}(∂²zᵃᵃᶠ_top, grid_gpu, ζ)
ζ_zz_face = Field(ζ_zz_op_face)
compute!(ζ_zz_face)

#=@inline δzᵃᵃᶠ(i, j, k, grid, u) = @inbounds u[i, j, k+1] - u[i, j, k]
ζ_zz_op_face = KernelFunctionOperation{Face, Face, Face}(∂²zᵃᵃᶜ, grid_gpu, ζ)
ζ_zz_face = Field(ζ_zz_op_face)
compute!(ζ_zz_face)=#

ζ_zz_new_cpu = Field{Face, Face, Center}(grid_cpu)
CUDA.@allowscalar ζ_zz_new_cpu .= (ζ_zz_new.data)
nothing

Δz = 10/Nz
zᶜ = 0.25:Δz:9.75
zᶠ = 0:Δz:10
CUDA.@allowscalar ζ_zz_arr = [ζ_zz_face[2,1,n] for n = 1 : length(zᶠ)]
CUDA.@allowscalar ζ_arr = [ζ[2,1,n] for n = 1 : length(zᶜ)]

ζᶜ_rep = [0; 0; ζ_arr[1]; ζ_arr; ζ_arr[end]; 0; 0]
ζ_zzᶠ_rep = [0; 0; 0; [(3ζᶜ_rep[k+1] - 7ζᶜ_rep[k] + 5ζᶜ_rep[k-1] - ζᶜ_rep[k-2]) / (2Δz^2) for k = 3 : length(zᶠ)+2]; 0; 0; 0]
ζ_zzᶜ_rep = [0; 0; 0; [(ζᶜ_rep[k+1] - 2ζᶜ_rep[k] + ζᶜ_rep[k-1]) / (Δz^2) for k = 4 : length(zᶜ)+3]; 0; 0; 0]
ζ_zzᶜ_rep[[3, end-2]] = ζ_zzᶜ_rep[[4, end-3]]

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, ζᶜ_rep[4:end-3], zᶜ)
#lines!(ax, ζ_zzᶠ_rep[6:end-3], zᶠ[3:end])
lines!(ax, ζ_zzᶜ_rep[4:end-3], zᶜ)
lines!(ax, ζ_zz_arr[3:end], zᶠ[3:end])
#lines!(ax, -k^2*f.(zᶜ), zᶜ)
lines!(ax, -k^2*f.(zᶠ), zᶠ)
display(fig)

# f_zz(0) = 1/2Δz² * [3f(Δz/2) - 7f(-Δz/2) + 5f(-3Δz/2) - f(-5Δz/2)] + 𝒪(Δz²)
# consider calculating ζ_zzᵃᵃᶠ

Δz = 10/Nz
fig = Figure()
ax = Axis(fig[1, 1])
zᶜ = -2.5Δz : Δz : 10+2.5Δz
CUDA.@allowscalar ζ_arr = [ζ[2,1,n] for n = -2 : length(zᶜ)-3]
∇ₕ²(𝑓) = ∂x(∂x(𝑓)) + ∂y(∂y(𝑓))
∇²δ = ∇ₕ²(δ)
CUDA.@allowscalar ∇²δ_arr = [∇²δ[5,2,n] for n = -2 : length(zᶜ)-3]
lines!(ax, ζ_arr, zᶜ)
lines!(ax, ∇²δ_arr, zᶜ)
lines!(ax, f.(zᶜ), zᶜ)
display(fig)

#=
fig = Figure()
ax = Axis(fig[1, 1])
z = -0.5Δz : Δz : 10+0.5Δz
lines!(ax, ζ[2,1,0:Nz+1], z)
lines!(ax, [(2ζ[2,1,0]-5ζ[2,1,1]+4ζ[2,1,2]-ζ[2,1,3])/Δz^2; [ζ_zz[2,1,n] for n = 1:Nz]; (-ζ[2,1,Nz-2]+4ζ[2,1,Nz-1]-5ζ[2,1,Nz]+2ζ[2,1,Nz+1])/Δz^2], z)
lines!(ax, -k^2*f.(z), z)
display(fig)=#

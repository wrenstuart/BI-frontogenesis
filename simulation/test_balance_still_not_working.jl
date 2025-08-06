# NOT RESOLVED

using Printf
using Oceananigans
using Oceananigans.TurbulenceClosures
using Oceananigans.Operators
using Oceananigans.BoundaryConditions: fill_halo_regions!
using CUDA
using StructArrays
using Oceananigans.Fields
using Oceananigans.Architectures: arch_array
using Unroll

using Oceananigans.Models.NonhydrostaticModels

label = "test_extra_visc_low_res"

include("inputs/" * label * ".jl")
include("../QOL.jl")
include("../instabilities/modes.jl")
include("tendies2.jl")

function physical_quantities_from_inputs(Ri, s)

    # Get the dimensional parameters of the problem
    p = get_scales(Ri, s)

    # Set the viscosities

    # Set the domain size
    Lx = 2 * 2π * p.L * 0.4^0.5 # Zonal extent, set to 2 wavelengths of the most unstable mode
    Ly = Lx                     # Meridional extent
    Lz = p.H                    # Vertical extent

    # Set relative amplitude for random velocity perturbation

    B₀(x, y, z, t) = p.M² * y + p.N² * z    # Buoyancy
    U₀(x, y, z, t) = -p.M²/p.f * (z + Lz)   # Zonal velocity

    # Set the initial perturbation conditions, a random velocity perturbation
    uᵢ, vᵢ, wᵢ, bᵢ = generate_ic(Ri, Lx, p.U, N = 10)

    u_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0),
                                    bottom = GradientBoundaryCondition(0.0))
    # This is applied to the perturbation u, not including the background U₀
    # I think this is not strictly necessary, as that is the default boundary
    # condition for the velocity field, but it is good to be explicit
    BCs = (u = u_bcs,)

    return p, (x = Lx, y = Ly, z = Lz), (u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ), (U = U₀, B = B₀), BCs

end

struct MyParticle

    x::Float64
    y::Float64
    z::Float64

    ζ::Float64
    δ::Float64
    ζ_tendency::Float64
    ζ_cor::Float64
    ζ_visc::Float64
    ζ_err::Float64
    F_ζ_hor::Float64
    F_ζ_vrt::Float64
    ζ_adv::Float64
    ζ_h_adv::Float64

end

@info label
dir = "raw_data/" * label
if isdir(dir)
    throw("Output directory for label " * label * " already exists")
else
    mkdir(dir)
    @info "Created director for simulation with label " * label
end

params = sim_params()
resolution = params.res

phys_params, domain, ic, background, BCs = physical_quantities_from_inputs(params.Ri, params.s)
f = phys_params.f

# Set the time-stepping parameters
max_Δt = 0.4 * pi / (phys_params.N²^0.5)
duration = 10 / real(least_stable_mode(params.Ri, 4π/domain.x, 0, rate_only = true))
if params.short_duration
    duration = duration / 20
end

# Build the grid
if params.GPU
    grid = RectilinearGrid(GPU(), size = resolution, x = (0, domain.x), y = (0, domain.y), z = (-domain.z, 0), topology = (Periodic, Periodic, Bounded))
else
    grid = RectilinearGrid(size = resolution, x = (0, domain.x), y = (0, domain.y), z = (-domain.z, 0), topology = (Periodic, Periodic, Bounded))
end

# Set the diffusivities and background fields
B_field = BackgroundField(background.B)
U_field = BackgroundField(background.U)
empty_field = BackgroundField(0.0)
if sim_params().horizontal_hyperviscosity
    diff_h = ScalarBiharmonicDiffusivity(Oceananigans.TurbulenceClosures.HorizontalFormulation(), Float64, ν = params.ν_h, κ = params.ν_h)
else
    diff_h = HorizontalScalarDiffusivity(ν = params.ν_h, κ = params.ν_h)
end
diff_v = VerticalScalarDiffusivity(ν = params.ν_v, κ = params.ν_v)

# Introduce Lagrangian particles in an n × n grid
n = 20
x₀ = Array{Float64, 1}(undef, n^2)
y₀ = Array{Float64, 1}(undef, n^2)
@unroll for i = 0 : n^2-1
    x₀[i+1] = domain.x * (i % n) / n
    y₀[i+1] = domain.y * (i ÷ n) / n
end
if params.GPU
    x₀, y₀ = CuArray.([x₀, y₀])
end
O = params.GPU ? () -> CuArray(zeros(n^2)) : () -> zeros(n^2)
particles = StructArray{MyParticle}((x₀, y₀, O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O()))

@inline function ∂²zᵃᵃᶠ_top(i, j, k, grid, u)
    δz² = Δzᶠᶠᶜ(i, j, k, grid)^2
    δ²u = (3u[i, j, k] - 7u[i, j, k-1] + 5u[i, j, k-2] - u[i, j, k-3]) / 2
    return δ²u/δz²
end
∇ₕ²(𝑓) = ∂x(∂x(𝑓)) + ∂y(∂y(𝑓))

# Extract fundamental variable fields:
velocities = VelocityFields(grid)
u, v, w = velocities
tracers = TracerFields((:b,), grid)
b, = tracers
pHY′ = CenterField(grid)
pNHS = CenterField(grid)

# Intermediary terms:
p = pHY′ + pNHS
u_y = ∂y(u)
v_y = ∂y(v)

# To store as auxiliary fields
u_x = ∂x(u)
v_x = ∂x(v)
M²_on_f = phys_params.M²/f
u_z = ∂z(u) - M²_on_f   # These five
v_z = ∂z(v)             # vanish or are
w_x = ∂x(w)             # constant at the
w_y = ∂y(w)             # boundaries, but not
b_z = ∂z(b)             # in the interior
b_x = ∂x(b)
b_y = ∂y(b) + phys_params.M²
ζ = v_x - u_y
ζ_z = ∂z(ζ)
δ = u_x + v_y
∇ₕ²ζ = ∇ₕ²(ζ)
∇ₕ²δ = ∇ₕ²(δ)
ζ_zz = ∂z(∂z(ζ))
δ_zz = ∂z(∂z(δ))
fζ_g = ∇ₕ²(p)
ζ_zz_op = KernelFunctionOperation{Face, Face, Face}(∂²zᵃᵃᶠ_top, grid, ζ)
δ_zz_op = KernelFunctionOperation{Center, Center, Face}(∂²zᵃᵃᶠ_top, grid, δ)
b_xzz_op = KernelFunctionOperation{Face, Center, Face}(∂²zᵃᵃᶠ_top, grid, b_x)
b_yzz_op = KernelFunctionOperation{Center, Face, Face}(∂²zᵃᵃᶠ_top, grid, b_y)
my_ζ_zz = Field(ζ_zz_op)
#my_δ_zz = Field(δ_zz_op)
b_xzz = Field(b_xzz_op)
b_yzz = Field(b_yzz_op)

my_ζ_visc = params.ν_h * ∇ₕ²ζ + params.ν_v * my_ζ_zz
my_ζ_cor = -ζ*δ
my_ζ_hor = -f*δ

u_back = Field{Face, Center, Center}(grid)
v_back = Field{Center, Face, Center}(grid)
w_back = Field{Center, Center, Face}(grid)
b_back = Field{Center, Center, Center}(grid)
set!(u_back, (x, y, z) -> background.U(x, y, z, 0))
set!(b_back, (x, y, z) -> background.B(x, y, z, 0))
fill_halo_regions!(u_back)
fill_halo_regions!(b_back)
background_fields = (velocities = (u = u_back, v = v_back, w = w_back),
                     tracers = (b_back))
closure = (diff_h, diff_v)
diffusivities = ((ν = params.ν_h, κ = params.ν_h), (ν = params.ν_v, κ = params.ν_v))

other_args = (advection_scheme = params.advection_scheme(),
              coriolis = FPlane(f = f),
              closure = closure,
              buoyancy = Buoyancy(model = BuoyancyTracer()),
              background_fields = background_fields,
              velocities = velocities,
              tracers = tracers,
              diffusivities = diffusivities,
              hydrostatic_pressure = pHY′)

@inline u_tendency_op = KernelFunctionOperation{Face, Center, Center}(u_tendency_func, grid, other_args)
@inline u_cor_op      = KernelFunctionOperation{Face, Center, Center}(u_cor_func,      grid, other_args)
@inline u_visc_op     = KernelFunctionOperation{Face, Center, Center}(u_visc_func,     grid, other_args)
@inline u_err_op      = KernelFunctionOperation{Face, Center, Center}(u_err_func,      grid, other_args)
@inline u_div𝐯_op     = KernelFunctionOperation{Face, Center, Center}(u_div𝐯_func,     grid, other_args)
@inline v_tendency_op = KernelFunctionOperation{Center, Face, Center}(v_tendency_func, grid, other_args)
@inline v_cor_op      = KernelFunctionOperation{Center, Face, Center}(v_cor_func,      grid, other_args)
@inline v_visc_op     = KernelFunctionOperation{Center, Face, Center}(v_visc_func,     grid, other_args)
@inline v_err_op      = KernelFunctionOperation{Center, Face, Center}(v_err_func,      grid, other_args)
@inline v_div𝐯_op     = KernelFunctionOperation{Center, Face, Center}(v_div𝐯_func,     grid, other_args)
@inline my_u_div𝐯_op  = KernelFunctionOperation{Face, Center, Center}(my_u_div𝐯_func,  grid, other_args)
@inline my_v_div𝐯_op  = KernelFunctionOperation{Center, Face, Center}(my_v_div𝐯_func,  grid, other_args)
u_tendency = Field(u_tendency_op)
u_cor      = Field(u_cor_op)
u_visc     = Field(u_visc_op)
u_err      = Field(u_err_op)
u_div𝐯     = Field(u_div𝐯_op)
v_tendency = Field(v_tendency_op)
v_cor      = Field(v_cor_op)
v_visc     = Field(v_visc_op)
v_err      = Field(v_err_op)
v_div𝐯     = Field(v_div𝐯_op)
ζ_tendency = ∂x(v_tendency) - ∂y(u_tendency)
ζ_cor      = ∂x(v_cor)      - ∂y(u_cor)
ζ_visc     = ∂x(v_visc)     - ∂y(u_visc)
ζ_err      = ∂x(v_err)      - ∂y(u_err)
# ζ_div𝐯     = ∂x(v_div𝐯)     - ∂y(u_div𝐯)        # 𝐳̂⋅∇×(∇⋅(𝐮𝐮))
compute!(ζ_tendency)
compute!(ζ_cor)
compute!(ζ_visc)
compute!(ζ_err)
my_u_div𝐯  = Field(my_u_div𝐯_op)
my_v_div𝐯  = Field(my_v_div𝐯_op)

#=@inline ζ_tendency_op = KernelFunctionOperation{Face, Face, Center}(ζ_tendency_func, grid, other_args)
ζ_tendency = Field(ζ_tendency_op)
@inline ζ_cor_op = KernelFunctionOperation{Face, Face, Center}(ζ_cor_func, grid, other_args)
ζ_cor = Field(ζ_cor_op)
@inline ζ_visc_op = KernelFunctionOperation{Face, Face, Center}(ζ_visc_func, grid, other_args)
ζ_visc = Field(ζ_visc_op)
@inline ζ_err_op = KernelFunctionOperation{Face, Face, Center}(ζ_err_func, grid, other_args)
ζ_err = Field(ζ_err_op)=#
@inline F_ζ_hor_op = KernelFunctionOperation{Face, Face, Center}(F_ζ_hor_func, grid, other_args)
@inline F_ζ_vrt_op = KernelFunctionOperation{Face, Face, Center}(F_ζ_vrt_func, grid, other_args)
@inline ζ_adv_op = KernelFunctionOperation{Face, Face, Center}(ζ_adv_func, grid, other_args)
@inline ζ_h_adv_op = KernelFunctionOperation{Face, Face, Center}(ζ_h_adv_func, grid, other_args)
F_ζ_hor = Field(F_ζ_hor_op)
F_ζ_vrt = Field(F_ζ_vrt_op)
ζ_adv = Field(ζ_adv_op)
ζ_h_adv = Field(ζ_h_adv_op)

auxiliary_fields = (; ζ, δ, ζ_tendency, ζ_cor, ζ_visc, ζ_err, F_ζ_hor, F_ζ_vrt, ζ_adv, ζ_h_adv)
drifter_fields = auxiliary_fields

function fix_particle_below_surface(lagrangian_particles, model, Δt)
    lagrangian_particles.properties.z .= -domain.z / 2resolution[3]
end
if params.fix_drifters_below_surface
    lagrangian_drifters = LagrangianParticles(particles; tracked_fields = drifter_fields, dynamics = fix_particle_below_surface)
else
    lagrangian_drifters = LagrangianParticles(particles; tracked_fields = drifter_fields)
end

# Remember to use CuArray instead of regular Array when storing particle locations and properties on the GPU?????

# Build the model
model = NonhydrostaticModel(; grid,
          advection = params.advection_scheme(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
          timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
          tracers = tracers,  # Set the name(s) of any tracers; here, b is buoyancy
          velocities = velocities,
          auxiliary_fields = auxiliary_fields,
          pressures = (; pHY′, pNHS),
          buoyancy = Buoyancy(model = BuoyancyTracer()), # this tells the model that b will act as the buoyancy (and influence momentum)
          background_fields = (b = B_field, u = U_field),
          coriolis = coriolis = FPlane(f = f),
          closure = (diff_h, diff_v),
          boundary_conditions = BCs,
          particles = lagrangian_drifters)

# Set initial conditions
set!(model, u = ic.u, v = ic.v, w = ic.w, b = ic.b)

# Build the simulation
simulation = Simulation(model, Δt = minimum([max_Δt/10, phys_params.T/100]), stop_time = duration)

# ### The `TimeStepWizard`
#
# The TimeStepWizard manages the time-step adaptively, keeping the
# Courant-Freidrichs-Lewy (CFL) number close to `1.0` while ensuring
# the time-step does not increase beyond the maximum allowable value
wizard = TimeStepWizard(cfl = 0.5, max_change = 1.1, max_Δt = max_Δt)

# Still some numerical noise at CFL 0.1 for Ri = 10⁴, but none for CFL = 0.05

# A "Callback" pauses the simulation after a specified number of timesteps and calls a function (here the timestep wizard to update the timestep)
# To update the timestep more or less often, change IterationInterval in the next line
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# ### A progress messenger
# We add a callback that prints out a helpful bprogress message while the simulation runs.

start_time = time_ns()

progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(1e-9 * (time_ns() - start_time)),
                        sim.Δt,
                        AdvectiveCFL(sim.Δt)(sim.model))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# ### Output

u = Field(model.velocities.u + model.background_fields.velocities.u)    # Unpack velocity `Field`s
v = Field(model.velocities.v)
w = Field(model.velocities.w)
b = Field(model.tracers.b + model.background_fields.tracers.b)          # Extract the buoyancy and add the background field
b_pert = Field(model.tracers.b)
p = Field(model.pressures.pNHS + model.pressures.pHY′)

# Compute y-averages 𝐮̅(x,z) and b̅(x,z)
u̅ = Field(Average(u, dims = 2))
v̅ = Field(Average(v, dims = 2))
w̅ = Field(Average(w, dims = 2))
b̅ = Field(Average(b, dims = 2))
ℬ = Field(w * b_pert)
avg_ℬ = Field(Average(ℬ, dims = 2))

# Output Lagrangian particles
filename = dir * "/particles"
simulation.output_writers[:particles] =
    JLD2OutputWriter(model, (particles = model.particles,),
                            filename = filename * ".jld2",
                            schedule = TimeInterval(phys_params.T/30),
                            overwrite_existing = true)

# Output the slice y = 0
#filename = dir * "/BI_xz"
#simulation.output_writers[:xz_slices] =
#    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, fζ_g),
#                            filename = filename * ".jld2",
#                            indices = (:, 1, :),
#                            schedule = TimeInterval(phys_params.T/100),
#                            overwrite_existing = true)

# Output the slice z = 0
filename = dir * "/BI_xy"
simulation.output_writers[:xy_slices] =
    JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, ζ_tendency, ζ_cor, ζ_visc, ζ_err, F_ζ_hor, F_ζ_vrt, ζ_adv, ζ_h_adv, u_div𝐯, v_div𝐯, my_u_div𝐯, my_v_div𝐯),
                            filename = filename * ".jld2",
                            indices = (:, :, resolution[3]),
                            schedule = TimeInterval(phys_params.T/30),
                            overwrite_existing = true)

nothing # hide

# Now, run the simulation
@printf("Simulation will last %s\n", prettytime(duration))
run!(simulation)
# Functions for running an individual simulation

using Oceananigans
using Printf
using Oceananigans.TurbulenceClosures
using CUDA

using StructArrays
using Oceananigans.Fields
using Oceananigans.Architectures: arch_array

include("../QOL.jl")
include("../instabilities/modes.jl")

function physical_quantities_from_inputs(Ri, s)

    # Get the dimensional parameters of the problem
    p = get_scales(Ri, s)

    # Set the viscosities

    # Set the domain size
    Lx = 2 * 2*pi * p.L * 0.4^0.5   # Zonal extent, set to 2 wavelengths of the most unstable mode
    Ly = Lx                         # Meridional extent
    Lz = p.H                        # Vertical extent

    # Set relative amplitude for random velocity perturbation
    kick = 0.05 * p.U

    B₀(x, y, z, t) = p.M² * y + p.N² * z    # Buoyancy
    U₀(x, y, z, t) = -p.M²/p.f * (z + Lz)   # Zonal velocity

    # Set the initial perturbation conditions, a random velocity perturbation
    uᵢ, vᵢ, wᵢ, bᵢ = generate_ic(Ri, Lx, p.U)

    u_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(0.0),
                                    bottom = GradientBoundaryCondition(0.0))
    # This is applied to the perturbation u, not including the background U₀
    BCs = (u = u_bcs,)

    return p, (x = Lx, y = Ly, z = Lz), (u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ), (U = U₀, B = B₀), BCs

end

struct MyParticle

    x::Float64
    y::Float64
    z::Float64

    u::Float64
    v::Float64
    w::Float64
    b::Float64
    p::Float64

    ζ::Float64
    δ::Float64
    u_x::Float64
    v_x::Float64
    u_z::Float64
    v_z::Float64
    w_x::Float64
    w_y::Float64
    b_x::Float64
    b_y::Float64
    b_z::Float64
    u_g::Float64
    v_g::Float64
    ζ_g::Float64    # Equivalently, ∇²Φ

    # For diffusive terms
    ζ_zz::Float64
    ∇ₕ²ζ::Float64
    δ_zz::Float64
    ∇ₕ²δ::Float64
    u_xzz::Float64
    ∇ₕ²u_x::Float64
    v_xzz::Float64
    ∇ₕ²v_x::Float64
    u_zzz::Float64
    ∇ₕ²u_z::Float64
    v_zzz::Float64
    ∇ₕ²v_z::Float64
    b_xzz::Float64
    ∇ₕ²b_x::Float64
    b_yzz::Float64
    ∇ₕ²b_y::Float64
    
    # Probably not necessary but just in case
    w_xzz::Float64
    ∇ₕ²w_x::Float64
    w_yzz::Float64
    ∇ₕ²w_y::Float64
    b_zzz::Float64
    ∇ₕ²b_z::Float64

end

function run_sim(params, label)

    resolution = params.res

    @info label

    phys_params, domain, ic, background, BCs = physical_quantities_from_inputs(params.Ri, params.s)
    f = phys_params.f

    # Set the time-stepping parameters
    max_Δt = 0.4 * pi / (phys_params.N²^0.5)
    duration = 1 / real(least_stable_mode(params.Ri, 4π/domain.x, 0, rate_only = true))

    # Build the grid
    if params.GPU
        grid = RectilinearGrid(GPU(), size = resolution, x = (0, domain.x), y = (0, domain.y), z = (-domain.z, 0), topology = (Periodic, Periodic, Bounded))
    else
        grid = RectilinearGrid(size = resolution, x = (0, domain.x), y = (0, domain.y), z = (-domain.z, 0), topology = (Periodic, Periodic, Bounded))
    end

    # Set the diffusivities and background fields
    B_field = BackgroundField(background.B)
    U_field = BackgroundField(background.U)
    if sim_params().horizontal_hyperviscosity
        diff_h = ScalarBiharmonicDiffusivity(Oceananigans.TurbulenceClosures.HorizontalFormulation(), Float64, ν = params.ν_h, κ = params.ν_h)
    else
        diff_h = HorizontalScalarDiffusivity(ν = params.ν_h, κ = params.ν_h)
    end
    diff_v = VerticalScalarDiffusivity(ν = params.ν_v, κ = params.ν_v)

    # Introduce Lagrangian particles in an n × n grid
    n = 10
    x₀ = [domain.x * (i % n) / n for i = 0 : n^2-1]
    y₀ = [domain.y * (i ÷ n) / n for i = 0: n^2-1]
    Os = zeros(n^2)
    if params.GPU
        x₀, y₀, Os = CuArray.([x₀, y₀, Os])
    end
    particles = StructArray{MyParticle}((x₀, y₀, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os, Os))

    tracers = TracerFields((:b,), grid)
    b = tracers[1]
    velocities = VelocityFields(grid)
    u, v, w = velocities
    pHY′ = CenterField(grid)
    pNHS = CenterField(grid)
    
    # To calculate (itermediary):
    p = pHY′ + pNHS
    u_y = ∂y(u)
    w_z = ∂z(w)

    # To store as an auxillary field:
    u_x = ∂x(u)
    v_x = ∂x(v)
    u_z = ∂z(u)
    v_z = ∂z(v)
    w_x = ∂x(w)
    w_y = ∂y(w)
    b_x = ∂x(b)
    b_y = ∂y(b)
    b_z = ∂z(b)
    ζ = v_x - u_y
    δ = -w_z
    u_g = -∂y(p)/f
    v_g = ∂x(p)/f
    ζ_g = ∂x(v_g) - ∂y(u_g)     # Equivalently, ∇²Φ

    ∇ₕ²(𝑓) = ∂x(∂x(𝑓)) + ∂y(∂y(𝑓))

    ζ_zz = ∂z(∂z(ζ))
    ∇ₕ²ζ = ∇ₕ²(ζ)
    δ_zz = ∂z(∂z(δ))
    ∇ₕ²δ = ∇ₕ²(δ)

    u_xzz = ∂z(∂z(u_x))
    ∇ₕ²u_x = ∇ₕ²(u_x)
    v_xzz = ∂z(∂z(v_x))
    ∇ₕ²v_x = ∇ₕ²(v_x)
    u_zzz = ∂z(∂z(u_z))
    ∇ₕ²u_z = ∇ₕ²(u_z)
    v_zzz = ∂z(∂z(v_z))
    ∇ₕ²v_z = ∇ₕ²(v_z)
    b_xzz = ∂z(∂z(b_x))
    ∇ₕ²b_x = ∇ₕ²(b_x)
    b_yzz = ∂z(∂z(b_y))
    ∇ₕ²b_y = ∇ₕ²(b_y)
    
    w_xzz = ∂z(∂z(w_x))
    ∇ₕ²w_x = ∇ₕ²(w_x)
    w_yzz = ∂z(∂z(w_y))
    ∇ₕ²w_y = ∇ₕ²(w_y)
    b_zzz = ∂z(∂z(b_z))
    ∇ₕ²b_z = ∇ₕ²(b_z)
    #######################Look into D(𝑁²)/D𝑡?

    #=∇ₕ²(𝑓) = ∂x(∂x(𝑓)) + ∂y(∂y(𝑓))
    F_hor_ζ = -δ * ζ
    F_cor_ζ = -f * δ
    F_ver_ζ = u_z * w_y - v_z * w_x
    V_mix_ζ = b#params.ν_v * ∂z(∂z(ζ))
    H_dif_ζ = b#params.ν_h * ∇ₕ²(ζ)
    F_hor_δ = -(u_x*u_x + 2u_y*v_x + v_y*v_y)
    F_hor_δ_approx = -δ^2
    F_cor_δ = f * ζ
    F_ver_δ = -(u_z*w_x + v_z*w_y)
    F_prs_δ = -f * ζ_g
    V_mix_δ = b#params.ν_v * ∂z(∂z(δ))
    H_dif_δ = b#params.ν_h * ∇ₕ²(δ)
    # Might not need these, but just in case:
    V_mix_b = b#params.ν_v * (b_x*∂z(∂z(b_x)) + b_y*∂z(∂z(b_y)))
    H_dif_b = b#params.ν_h * (b_x*∇ₕ²(b_x) + b_y*∇ₕ²(b_y))
    V_mix_u = b#params.ν_v * (u_x*∂z(∂z(u_x)) + u_y*∂z(∂z(u_y)) + v_x*∂z(∂z(v_x)) + v_y*∂z(∂z(v_y)))
    H_dif_u = b#params.ν_h * (u_x*∇ₕ²(u_x) + u_y*∇ₕ²(u_y) + v_x*∇ₕ²(v_x) + v_y*∇ₕ²(v_y))=#

    tracked_fields = (; u, v, w, b, p, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, u_g, v_g, ζ_g, ζ_zz, ∇ₕ²ζ, δ_zz, ∇ₕ²δ, u_xzz, ∇ₕ²u_x, v_xzz, ∇ₕ²v_x, u_zzz, ∇ₕ²u_z, v_zzz, ∇ₕ²v_z, b_xzz, ∇ₕ²b_x, b_yzz, ∇ₕ²b_y, w_xzz, ∇ₕ²w_x, w_yzz, ∇ₕ²w_y, b_zzz, ∇ₕ²b_z)
    lagrangian_drifters = LagrangianParticles(particles; tracked_fields = tracked_fields)          

    # "Remember to use CuArray instead of regular Array when storing particle locations and properties on the GPU"?????

    # Build the model
    model = NonhydrostaticModel(; grid,
              advection = params.advection_scheme(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
              timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
              tracers = (; b),  # Set the name(s) of any tracers; here, b is buoyancy
              velocities = velocities,
              auxiliary_fields = (; p, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, u_g, v_g, ζ_g, ζ_zz, ∇ₕ²ζ, δ_zz, ∇ₕ²δ, u_xzz, ∇ₕ²u_x, v_xzz, ∇ₕ²v_x, u_zzz, ∇ₕ²u_z, v_zzz, ∇ₕ²v_z, b_xzz, ∇ₕ²b_x, b_yzz, ∇ₕ²b_y, w_xzz, ∇ₕ²w_x, w_yzz, ∇ₕ²w_y, b_zzz, ∇ₕ²b_z, ),
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
    # We add a callback that prints out a helpful progress message while the simulation runs.

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
    p = Field(model.auxiliary_fields.p)    # (should i do + .pHY′?) # (ignoring the background pressure field, M²y(z+H) + N²z²/2)
    ζ = Field(model.auxiliary_fields.ζ)
    δ = Field(model.auxiliary_fields.δ)    # The horizontal divergence
    u_x = Field(model.auxiliary_fields.u_x)
    v_x = Field(model.auxiliary_fields.v_x)
    u_z = Field(model.auxiliary_fields.u_z)
    v_z = Field(model.auxiliary_fields.v_z)
    w_x = Field(model.auxiliary_fields.w_x)
    w_y = Field(model.auxiliary_fields.w_y)
    u_g = Field(model.auxiliary_fields.u_g)
    v_g = Field(model.auxiliary_fields.v_g)
    ζ_g = Field(model.auxiliary_fields.ζ_g)
    b_x = Field(model.auxiliary_fields.b_x)
    b_y = Field(model.auxiliary_fields.b_y)
    b_z = Field(model.auxiliary_fields.b_z)

    # Compute y-averages 𝐮̅(x,z) and b̅(x,z)
    u̅ = Field(Average(u, dims = 2))
    v̅ = Field(Average(v, dims = 2))
    w̅ = Field(Average(w, dims = 2))
    b̅ = Field(Average(b, dims = 2))
    ℬ = Field(w * b_pert)
    avg_ℬ = Field(Average(ℬ, dims = 2))

    # Output Lagrangian particles
    filename = "raw_data/" * label * "_particles"
    simulation.output_writers[:particles] =
        JLD2OutputWriter(model, (particles = model.particles,),
                                filename = filename * ".jld2",
                                schedule = TimeInterval(phys_params.T/100),
                                overwrite_existing = true)

    # Output the slice y = 0
    filename = "raw_data/" * label * "_BI_xz"
    simulation.output_writers[:xz_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, u_g, v_g, ζ_g),
                                filename = filename * ".jld2",
                                indices = (:, 1, :),
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    # Output the slice z = 0
    filename = "raw_data/" * label * "_BI_xy"
    simulation.output_writers[:xy_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, u_g, v_g, ζ_g),
                                filename = filename * ".jld2",
                                indices = (:, :, resolution[3]),
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    # Output the slice x = 0
    filename = "raw_data/" * label * "_BI_yz"
    simulation.output_writers[:yz_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, u_g, v_g, ζ_g),
                                filename = filename * ".jld2",
                                indices = (1, :, :),
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    # Output a horizontal slice in the middle (verticall speaking)
    filename = "raw_data/" * label * "_BI_xy_mid"
    simulation.output_writers[:xy_slices_mid] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, u_g, v_g, ζ_g),
                                filename = filename * ".jld2",
                                indices = (:, :, Int64(round((resolution[3]+1) / 2))),
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    filename = "raw_data/" * label * "_BI_y-avg"
    simulation.output_writers[:xy_slices_mid] =
        JLD2OutputWriter(model, (; u̅, v̅, w̅, b̅, avg_ℬ),
                                filename = filename * ".jld2",
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    #=filename = "raw_data/" * label * "_full"
    simulation.output_writers[:full] =
        JLD2OutputWriter(model, (; u, v, w, b),
                                filename = filename * ".jld2",
                                schedule = TimeInterval(duration / 4),
                                overwrite_existing = true)=#

    nothing # hide

    # Now, run the simulation
    @printf("Simulation will last %s\n", prettytime(duration))
    run!(simulation)

end
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
    b::Float64
    ζ::Float64
    δ::Float64
end

function run_sim(params, label)

    resolution = params.res

    @info label

    p, domain, ic, background, BCs = physical_quantities_from_inputs(params.Ri, params.s)
    
    # Here, p holds the physical parameters

    # Set the time-stepping parameters
    max_Δt = 0.4 * pi / (p.N²^0.5)
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
    n = 5
    x₀ = [domain.x * (i % n) / n for i = 0 : n^2-1]
    y₀ = [domain.y * (i ÷ n) / n for i = 0: n^2-1]
    z₀ = zeros(n^2)
    b = z₀
    ζ = z₀
    δ = z₀
    if params.GPU
        x₀, y₀, z₀, b, ζ, δ = CuArray.([x₀, y₀, z₀, b, ζ, δ])
    end

    #=b = arch_array(CPU(), zeros(n^2))
    ζ = arch_array(CPU(), zeros(n^2))
    δ = arch_array(CPU(), zeros(n^2))=#
    particles = StructArray{MyParticle}((x₀, y₀, z₀, b, ζ, δ))

    #=function drifter_dynamics(particles, model, Δt)
        u = Field(model.velocities.u + model.background_fields.velocities.u)    # Unpack velocity `Field`s
        v = Field(model.velocities.v)
        ζ = Field(∂x(v) - ∂y(u))
        model.particles.properties.ζ .= ζ
    end=#

    function drifter_dynamics!(particles, model, Δt)
        u = Field(model.velocities.u + model.background_fields.velocities.u)    # Unpack velocity `Field`s
        v = Field(model.velocities.v)
        ζ = Field(∂x(v) - ∂y(u))
        N = length(particles.properties.x)
        x = particles.properties.x
        y = particles.properties.y
        z = particles.properties.z
        𝐱 = [(x[i], y[i], z[i]) for i = 1 : N]
        model.particles.properties.ζ .= [interpolate(𝐱[i], ζ) for i = 1 : N]
    end

    tracers = TracerFields((:b,), grid)
    b = tracers[1]

    # "Remember to use CuArray instead of regular Array when storing particle locations and properties on the GPU"?????

    # Build the model
    model = NonhydrostaticModel(; grid,
              advection = params.advection_scheme(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
              timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
              tracers = :b,  # Set the name(s) of any tracers, here b is buoyancy and c is a passive tracer (e.g. dye)
              buoyancy = Buoyancy(model = BuoyancyTracer()), # this tells the model that b will act as the buoyancy (and influence momentum)
              background_fields = (b = B_field, u = U_field),
              coriolis = coriolis = FPlane(f = p.f),
              closure = (diff_h, diff_v),
              boundary_conditions = BCs,
              )

    tracked_fields = model.tracers
    lagrangian_drifters = LagrangianParticles(particles; tracked_fields = (; b), dynamics = drifter_dynamics!)          

    model = NonhydrostaticModel(; grid,
              advection = params.advection_scheme(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
              timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
              tracers = :b,  # Set the name(s) of any tracers, here b is buoyancy and c is a passive tracer (e.g. dye)
              buoyancy = Buoyancy(model = BuoyancyTracer()), # this tells the model that b will act as the buoyancy (and influence momentum)
              background_fields = (b = B_field, u = U_field),
              coriolis = coriolis = FPlane(f = p.f),
              closure = (diff_h, diff_v),
              boundary_conditions = BCs,
              particles = lagrangian_drifters
              )

    # Set initial conditions
    set!(model, u = ic.u, v = ic.v, w = ic.w, b = ic.b)
    
    # Build the simulation
    simulation = Simulation(model, Δt = minimum([max_Δt/10, p.T/100]), stop_time = duration)

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
    p̃ = Field(model.pressures.pNHS)    # (should i do + .pHY′?) # (ignoring the background pressure field, M²y(z+H) + N²z²/2)

    # Now calculate the derivatives of 𝐮
    # Only 8 are needed, since ∇⋅𝐮 = 0
    ζ = Field(∂x(v) - ∂y(u))
    δ = Field(∂x(u) + ∂y(v))    # The horizontal divergence
    u_x = Field(∂x(u))
    v_x = Field(∂x(v))
    u_z = Field(∂z(u))
    v_z = Field(∂z(v))
    w_x = Field(∂x(w))
    w_y = Field(∂y(w))
    u_g = Field(∂y(p̃)/p.f + model.background_fields.velocities.u)
    v_g = Field(-∂x(p̃)/p.f)
    ζ_g = Field(∂x(v_g) - ∂y(u_g))

    # Also calculate derivatives of b
    b_x = Field(∂x(b))
    b_y = Field(Field(∂y(b_pert)) + p.M²)
    b_z = Field(∂z(b))

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
                                schedule = TimeInterval(p.T/20),
                                overwrite_existing = true)

    # Output the slice y = 0
    filename = "raw_data/" * label * "_BI_xz"
    simulation.output_writers[:xz_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, u_g, v_g, ζ_g),
                                filename = filename * ".jld2",
                                indices = (:, 1, :),
                                schedule = TimeInterval(p.T/20),
                                overwrite_existing = true)

    # Output the slice z = 0
    filename = "raw_data/" * label * "_BI_xy"
    simulation.output_writers[:xy_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, u_g, v_g, ζ_g),
                                filename = filename * ".jld2",
                                indices = (:, :, resolution[3]),
                                schedule = TimeInterval(p.T/20),
                                overwrite_existing = true)

    # Output the slice x = 0
    filename = "raw_data/" * label * "_BI_yz"
    simulation.output_writers[:yz_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, u_g, v_g, ζ_g),
                                filename = filename * ".jld2",
                                indices = (1, :, :),
                                schedule = TimeInterval(p.T/20),
                                overwrite_existing = true)

    # Output a horizontal slice in the middle (verticall speaking)
    filename = "raw_data/" * label * "_BI_xy_mid"
    simulation.output_writers[:xy_slices_mid] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, u_g, v_g, ζ_g),
                                filename = filename * ".jld2",
                                indices = (:, :, Int64(round((resolution[3]+1) / 2))),
                                schedule = TimeInterval(p.T/20),
                                overwrite_existing = true)

    filename = "raw_data/" * label * "_BI_y-avg"
    simulation.output_writers[:xy_slices_mid] =
        JLD2OutputWriter(model, (; u̅, v̅, w̅, b̅, avg_ℬ),
                                filename = filename * ".jld2",
                                schedule = TimeInterval(p.T/20),
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
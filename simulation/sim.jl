# Functions for running an individual simulation

using Oceananigans
using Printf

include("../QOL.jl")

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

    # Define the background fields
    #=B₀(x, y, z, t) = p.M² * y + p.N² * z    # Buoyancy
    U₀(x, y, z, t) = -p.M²/p.f * (z + Lz)   # Zonal velocity

    # Set the initial perturbation conditions, a random velocity perturbation
    uᵢ(x, y, z) = kick * randn()
    vᵢ(x, y, z) = kick * randn() 
    wᵢ(x, y, z) = kick * randn()
    bᵢ(x, y, z) = 0=#

    B₀(x, y, z, t) = p.M² * y   # Buoyancy
    U₀(x, y, z, t) = 0          # Zonal velocity

    # Set the initial perturbation conditions, a random velocity perturbation
    uᵢ(x, y, z) = kick * randn()
    vᵢ(x, y, z) = kick * randn() 
    wᵢ(x, y, z) = kick * randn()
    bᵢ(x, y, z) = p.N² * z

    return p, (x = Lx, y = Ly, z = Lz), (u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ), (U = U₀, B = B₀)

end

function run_sim(params)

    resolution = params.res
    label = params.label

    @info label

    p, domain, ic, background = physical_quantities_from_inputs(params.Ri, params.s)
    # Here, p holds the physical parameters

    # Set the time-stepping parameters
    max_Δt = minimum([p.T/10, 0.85 * pi / (p.N²^0.5)])
    duration = p.T * 40

    # Build the grid
    if params.GPU
        grid = RectilinearGrid(GPU(), size = resolution, x = (0, domain.x), y = (0, domain.y), z = (-domain.z, 0), topology = (Periodic, Periodic, Bounded))
    else
        grid = RectilinearGrid(size = resolution, x = (0, domain.x), y = (0, domain.y), z = (-domain.z, 0), topology = (Periodic, Periodic, Bounded))
    end

    # Set the diffusivities and background fields
    B_field = BackgroundField(background.B)
    U_field = BackgroundField(background.U)
    diff_h = HorizontalScalarDiffusivity(ν = params.ν_h, κ = params.ν_h)
    diff_v = HorizontalScalarDiffusivity(ν = params.ν_v, κ = params.ν_v)

    # Build the model
    model = NonhydrostaticModel(; grid,
              advection = params.advection_scheme(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
              timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
              tracers = :b,  # Set the name(s) of any tracers, here b is buoyancy and c is a passive tracer (e.g. dye)
              buoyancy = Buoyancy(model = BuoyancyTracer()), # this tells the model that b will act as the buoyancy (and influence momentum)
              background_fields = (b = B_field, u = U_field),
              coriolis = coriolis = FPlane(f = p.f),
              closure = (diff_h, diff_v)
              )
    
    # Set initial conditions
    set!(model, u = ic.u, v = ic.v, w = ic.w, b = ic.b)
    
    # Build the simulation
    simulation = Simulation(model, Δt = max_Δt/10, stop_time = duration)

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
    b_pert = Field(model.tracers.b)# - BackgroundField((x, y, z, t) -> p.N² * z))

    # Now calculate the derivatives of 𝐮
    # Only 8 are needed, since ∇⋅𝐮 = 0
    ζ₁ = Field(∂y(w) - ∂z(v))
    ζ₂ = Field(∂z(u) - ∂x(w))
    ζ₃ = Field(∂x(v) - ∂y(u))
    δ = Field(∂x(u) + ∂y(v))    # The horizontal divergence
    u_x = Field(∂x(u))
    v_x = Field(∂x(v))
    u_z = Field(∂z(u))
    v_z = Field(∂z(v))

    # Also calculate derivatives of b
    b_x = Field(∂x(b))
    b_y = Field(∂y(b))
    b_z = Field(∂z(b))

    # Compute y-averages 𝐮̅(x,z) and b̅(x,z)
    u̅ = Field(Average(u, dims = 2))
    v̅ = Field(Average(v, dims = 2))
    w̅ = Field(Average(w, dims = 2))
    b̅ = Field(Average(b, dims = 2))
    ℬ = Field(w * b_pert)
    avg_ℬ = Field(Average(ℬ, dims = 2))

    # Output the slice y = 0
    filename = "raw_data/" * label * "_BI_xz"
    simulation.output_writers[:xz_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ₁, ζ₂, ζ₃, δ, u_x, v_x, u_z, v_z, b_x, b_y, b_z),
                                filename = filename * ".jld2",
                                indices = (:, 1, :),
                                schedule = TimeInterval(p.T/40),
                                overwrite_existing = true)

    # Output the slice z = 0
    filename = "raw_data/" * label * "_BI_xy"
    simulation.output_writers[:xy_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ₁, ζ₂, ζ₃, δ, u_x, v_x, u_z, v_z, b_x, b_y, b_z),
                                filename = filename * ".jld2",
                                indices = (:, :, resolution[3]),
                                schedule = TimeInterval(p.T/40),
                                overwrite_existing = true)

    # Output the slice x = 0
    filename = "raw_data/" * label * "_BI_yz"
    simulation.output_writers[:yz_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ₁, ζ₂, ζ₃, δ, u_x, v_x, u_z, v_z, b_x, b_y, b_z),
                                filename = filename * ".jld2",
                                indices = (1, :, :),
                                schedule = TimeInterval(p.T/40),
                                overwrite_existing = true)

    # Output a horizontal slice in the middle (verticall speaking)
    filename = "raw_data/" * label * "_BI_xy_mid"
    simulation.output_writers[:xy_slices_mid] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ₁, ζ₂, ζ₃, δ, u_x, v_x, u_z, v_z, b_x, b_y, b_z),
                                filename = filename * ".jld2",
                                indices = (:, :, Int64(round((resolution[3]+1) / 2))),
                                schedule = TimeInterval(p.T/40),
                                overwrite_existing = true)

    filename = "raw_data/" * label * "_BI_y-avg"
    simulation.output_writers[:xy_slices_mid] =
        JLD2OutputWriter(model, (; u̅, v̅, w̅, b̅, avg_ℬ),
                                filename = filename * ".jld2",
                                schedule = TimeInterval(p.T/40),
                                overwrite_existing = true)

    filename = "raw_data/" * label * "_full"
    simulation.output_writers[:full] =
        JLD2OutputWriter(model, (; u, v, w, b),
                                filename = filename * ".jld2",
                                schedule = TimeInterval(duration / 4),
                                overwrite_existing = true)

    nothing # hide

    # Now, run the simulation
    @printf("Simulation will last %s\n", prettytime(duration))
    run!(simulation)

end
# Functions for running an individual simulation

using Printf
using Oceananigans
using Oceananigans.TurbulenceClosures
using Oceananigans.Operators
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

end

function run_sim(params, label)

    resolution = params.res

    @info label

    phys_params, domain, ic, background, BCs = physical_quantities_from_inputs(params.Ri, params.s)
    f = phys_params.f

    # Set the time-stepping parameters
    max_Δt = 0.4 * pi / (phys_params.N²^0.5)
    duration = 20 / real(least_stable_mode(params.Ri, 4π/domain.x, 0, rate_only = true))
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
    if sim_params().horizontal_hyperviscosity
        diff_h = ScalarBiharmonicDiffusivity(Oceananigans.TurbulenceClosures.HorizontalFormulation(), Float64, ν = params.ν_h, κ = params.ν_h)
    else
        diff_h = HorizontalScalarDiffusivity(ν = params.ν_h, κ = params.ν_h)
    end
    diff_v = VerticalScalarDiffusivity(ν = params.ν_v, κ = params.ν_v)

    # Introduce Lagrangian particles in an n × n grid
    n = 20
    x₀ = [domain.x * (i % n) / n for i = 0 : n^2-1]
    y₀ = [domain.y * (i ÷ n) / n for i = 0: n^2-1]
    if params.GPU
        x₀, y₀ = CuArray.([x₀, y₀])
    end
    O = params.GPU ? () -> CuArray(zeros(n^2)) : () -> zeros(n^2)
    particles = StructArray{MyParticle}((x₀, y₀, O()))#, O(), O()))#, O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O()))

    #=@inline δ²zᵃᵃᶠ_top(i, j, k, grid, u) = @inbounds (3u[i, j, k] - 7u[i, j, k-1] + 5u[i, j, k-2] - u[i, j, k-3]) / 2
    @inline ∂²zᵃᵃᶠ_top(i, j, k, grid, u) = @inbounds δ²zᵃᵃᶠ(i, j, k, grid, u) / Δzᶠᶠᶜ(i, j, k, grid)^2=#

    #=@inline function ∂²zᵃᵃᶠ_top(i, j, k, grid, u)

        Δz² = Δzᶠᶠᶜ(i, j, k, grid)^2
        δ²z = (3u[i, j, k] - 7u[i, j, k-1] + 5u[i, j, k-2] - u[i, j, k-3]) / 2
        return δ²z / Δz²
        
    end=#

    # Extract fundamental variable fields:
    velocities = VelocityFields(grid)
    u, v, w = velocities
    tracers = TracerFields((:b,), grid)
    b, = tracers
    pHY′ = CenterField(grid)
    pNHS = CenterField(grid)
    
    # Intermediary terms:
    p = pHY′ + pNHS
    
    lagrangian_drifters = LagrangianParticles(particles)

    # Remember to use CuArray instead of regular Array when storing particle locations and properties on the GPU?????

    # Build the model
    model = NonhydrostaticModel(; grid,
              advection = params.advection_scheme(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
              timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
              tracers = tracers,  # Set the name(s) of any tracers; here, b is buoyancy
              velocities = velocities,
              #auxiliary_fields = tracked_fields,
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
    p = Field(model.pressures.pNHS + model.pressures.pHY′)

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
        JLD2OutputWriter(model, (; u, v, w, b),
                                filename = filename * ".jld2",
                                indices = (:, 1, :),
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    # Output the slice z = 0
    filename = "raw_data/" * label * "_BI_xy"
    simulation.output_writers[:xy_slices] =
        JLD2OutputWriter(model, (; u, v, w, b),
                                filename = filename * ".jld2",
                                indices = (:, :, resolution[3]),
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    # Output the slice x = 0
    filename = "raw_data/" * label * "_BI_yz"
    simulation.output_writers[:yz_slices] =
        JLD2OutputWriter(model, (; u, v, w, b),
                                filename = filename * ".jld2",
                                indices = (1, :, :),
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    # Output a horizontal slice in the middle (verticall speaking)
    filename = "raw_data/" * label * "_BI_xy_mid"
    simulation.output_writers[:xy_slices_mid] =
        JLD2OutputWriter(model, (; u, v, w, b),
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
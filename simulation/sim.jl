# Functions for running an individual simulation

using Dates
using Printf
using Unroll
using CUDA
using StructArrays
using Oceananigans
using Oceananigans.TurbulenceClosures
using Oceananigans.Operators
using Oceananigans.Fields
using Oceananigans.Models.NonhydrostaticModels
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Architectures: arch_array

include("../QOL.jl")
include("../instabilities/modes.jl")
include("../io.jl")
include("tendies.jl")

function physical_quantities_from_inputs(Ri, s)

    # Get the dimensional parameters of the problem
    p = get_scales(Ri, s)

    # Set the domain size
    Lx = 2π * p.L * 0.4^0.5 # Zonal extent, set to 2 wavelengths of the most unstable mode
    Ly = Lx                     # Meridional extent
    Lz = p.H                    # Vertical extent

    B₀(x, y, z, t) = p.M² * y + p.N² * z        # Buoyancy
    U₀(x, y, z, t) = -p.M²/p.f * (z + Lz)       # Zonal velocity
    uᵢ, vᵢ, wᵢ, bᵢ = generate_ic(Ri, Lx, p.U)   # Initial conditions

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

    # Lagrangian ζ LHS:
    ζ_t::Float64            # 𝜁ₜ
    ζ_adv::Float64          # 𝐮⋅∇𝜁
    ζ_err::Float64          # 𝐳̂⋅∇×(∇⋅(𝐮𝐮)-𝐮⋅∇𝐮)
    # ^ vanishes in the continuum limit but non-0 when discretised
    # Lagrangian ζ RHS:
    F_ζ_hor::Float64        # -𝛿𝜁
    F_ζ_vrt::Float64        # -𝐳̂⋅(∇𝑤×𝚲)
    F_ζ_cor::Float64        # -𝛿𝑓
    ζ_h_visc::Float64       # νₕ∇ₕ²𝜁
    ζ_v_visc::Float64       # νᵥ∂²𝜁/∂𝑧²

    # Lagrangian δ LHS:
    δ_t::Float64            # 𝛿ₜ
    δ_adv::Float64          # 𝐮⋅∇𝛿
    δ_err::Float64          # ∇ₕ⋅∇(∇⋅(𝐮𝐮)-𝐮⋅∇𝐮) (c.f. ζ_err)
    # Lagrangian δ RHS:
    F_δ_hor::Float64        # -(∇ₕ𝐮ₕ):(∇ₕ𝐮ₕ)ᵀ
    F_δ_vrt::Float64        # -∇ₕ𝑤⋅𝚲
    F_δ_cor::Float64        # 𝑓𝜁
    F_δ_prs::Float64        # -∇ₕ²𝑝 ( = -𝑓𝜁_g)
    δ_h_visc::Float64       # νₕ∇ₕ²𝛿
    δ_v_visc::Float64       # νᵥ∂²𝛿/∂𝑧²

end

function run_sim(params, label; pickup = false)

    start_time = time_ns()          # (in nanoseconds)
    wall_time_limit = 5*60*60       # (in seconds)

    doubleoutput("Calculating physical parameters, including initial conditions", label)
    resolution = params.res
    phys_params, domain, ic, background, BCs = physical_quantities_from_inputs(params.Ri, params.s)
    f = phys_params.f

    doubleoutput("Calculating least stable mode", label)
    max_Δt = 0.4 * pi / (phys_params.N²^0.5)
    duration = 8 / real(least_stable_mode(params.Ri, 2π/domain.x, 0, rate_only = true))
    if params.short_duration
        duration = duration / 20
    end

    doubleoutput("Building the grid", label)
    if params.GPU
        grid = RectilinearGrid(GPU(), size = resolution, x = (0, domain.x), y = (0, domain.y), z = (-domain.z, 0), topology = (Periodic, Periodic, Bounded))
    else
        grid = RectilinearGrid(size = resolution, x = (0, domain.x), y = (0, domain.y), z = (-domain.z, 0), topology = (Periodic, Periodic, Bounded))
    end

    doubleoutput("Setting diffusivities and background fields", label)
    B_field = BackgroundField(background.B)
    U_field = BackgroundField(background.U)
    if sim_params().horizontal_hyperviscosity
        diff_h = ScalarBiharmonicDiffusivity(Oceananigans.TurbulenceClosures.HorizontalFormulation(), Float64, ν = params.ν_h, κ = params.ν_h)
    else
        diff_h = HorizontalScalarDiffusivity(ν = params.ν_h, κ = params.ν_h)
    end
    diff_v = VerticalScalarDiffusivity(ν = params.ν_v, κ = params.ν_v)

    doubleoutput("Generating particle array", label)
    x₀ = Array{Float64, 1}(undef, n_d^2)
    y₀ = Array{Float64, 1}(undef, n_d^2)
    @unroll for i = 0 : n_d^2-1
        x₀[i+1] = domain.x * (i % n_d) / n_d
        y₀[i+1] = domain.y * (i ÷ n_d) / n_d
    end
    if params.GPU
        x₀, y₀ = CuArray.([x₀, y₀])
    end
    O = params.GPU ? () -> CuArray(zeros(n_d^2)) : () -> zeros(n_d^2)
    particles = StructArray{MyParticle}((x₀, y₀, O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O()))

    doubleoutput("Building fundamental fields", label)
    u_back = Field{Face, Center, Center}(grid)
    v_back = Field{Center, Face, Center}(grid)
    w_back = Field{Center, Center, Face}(grid)
    b_back = Field{Center, Center, Center}(grid)
    set!(u_back, (x, y, z) -> background.U(x, y, z, 0))
    set!(b_back, (x, y, z) -> background.B(x, y, z, 0))
    fill_halo_regions!(u_back)
    fill_halo_regions!(b_back)
    velocities = VelocityFields(grid)
    tracers = TracerFields((:b,), grid)
    u_pert, v, w = velocities
    b_pert, = tracers
    u = u_pert + u_back
    b = b_pert + b_back
    pHY′ = CenterField(grid)
    pNHS = CenterField(grid)
    p = pHY′ + pNHS
    ζ = ∂x(v) - ∂y(u_pert)
    δ = ∂x(u_pert) + ∂y(v)

    background_fields = (velocities = (u = u_back, v = v_back, w = w_back),
                        tracers = (b_back))
    closure = (diff_h, diff_v)
    diffusivities = ((ν = params.ν_h, κ = params.ν_h), (ν = params.ν_v, κ = params.ν_v))
    other_args = (advection_scheme = params.advection_scheme(),
                coriolis = FPlane(f = f),
                closure = closure,
                buoyancy = BuoyancyTracer(),
                background_fields = background_fields,
                velocities = velocities,
                tracers = tracers,
                diffusivities = diffusivities,
                hydrostatic_pressure = pHY′,
                nonhydrostatic_pressure = pNHS)

    doubleoutput("Generating operators from kernels", label)
    @inline ζ_t      = Field(KernelFunctionOperation{Face, Face, Center}(ζ_t_func,      grid, other_args))
    @inline ζ_adv    = Field(KernelFunctionOperation{Face, Face, Center}(ζ_adv_func,    grid, other_args))
    @inline ζ_err    = Field(KernelFunctionOperation{Face, Face, Center}(ζ_err_func,    grid, other_args))
    @inline F_ζ_hor  = Field(KernelFunctionOperation{Face, Face, Center}(F_ζ_hor_func,  grid, other_args))
    @inline F_ζ_vrt  = Field(KernelFunctionOperation{Face, Face, Center}(F_ζ_vrt_func,  grid, other_args))
    @inline F_ζ_cor  = Field(KernelFunctionOperation{Face, Face, Center}(F_ζ_cor_func,  grid, other_args))
    @inline ζ_h_visc = Field(KernelFunctionOperation{Face, Face, Center}(ζ_h_visc_func, grid, other_args))
    @inline ζ_v_visc = Field(KernelFunctionOperation{Face, Face, Center}(ζ_v_visc_func, grid, other_args))
    @inline δ_t      = Field(KernelFunctionOperation{Center, Center, Center}(δ_t_func,      grid, other_args))
    @inline δ_adv    = Field(KernelFunctionOperation{Center, Center, Center}(δ_adv_func,    grid, other_args))
    @inline δ_err    = Field(KernelFunctionOperation{Center, Center, Center}(δ_err_func,    grid, other_args))
    @inline F_δ_hor  = Field(KernelFunctionOperation{Center, Center, Center}(F_δ_hor_func,  grid, other_args))
    @inline F_δ_vrt  = Field(KernelFunctionOperation{Center, Center, Center}(F_δ_vrt_func,  grid, other_args))
    @inline F_δ_cor  = Field(KernelFunctionOperation{Center, Center, Center}(F_δ_cor_func,  grid, other_args))
    @inline F_δ_prs  = Field(KernelFunctionOperation{Center, Center, Center}(F_δ_prs_func,  grid, other_args))
    @inline δ_h_visc = Field(KernelFunctionOperation{Center, Center, Center}(δ_h_visc_func, grid, other_args))
    @inline δ_v_visc = Field(KernelFunctionOperation{Center, Center, Center}(δ_v_visc_func, grid, other_args))
    fζ_g = - F_δ_prs

    drifter_fields = (; ζ, δ, ζ_t, ζ_adv, ζ_err, F_ζ_hor, F_ζ_vrt, F_ζ_cor, ζ_h_visc, ζ_v_visc, δ_t, δ_adv, δ_err, F_δ_hor, F_δ_vrt, F_δ_cor, F_δ_prs, δ_h_visc, δ_v_visc)

    doubleoutput("Generating particles", label)
    lagrangian_drifters = LagrangianParticles(particles; tracked_fields = drifter_fields)

    # Remember to use CuArray instead of regular Array when storing particle locations and properties on the GPU?????

    doubleoutput("Building model", label)
    model = NonhydrostaticModel(grid;
              advection = params.advection_scheme(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
              timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
              tracers = tracers,  # Set the name(s) of any tracers; here, b is buoyancy
              velocities = velocities,
              nonhydrostatic_pressure = pNHS,
              hydrostatic_pressure_anomaly = pHY′,
              buoyancy = BuoyancyTracer(), # this tells the model that b will act as the buoyancy (and influence momentum)
              background_fields = (b = B_field, u = U_field),
              coriolis = coriolis = FPlane(f = f),
              closure = (diff_h, diff_v),
              boundary_conditions = BCs,
              particles = lagrangian_drifters)

    if !pickup
        doubleoutput("Setting initial conditions", label)
        set!(model, u = ic.u, v = ic.v, w = ic.w, b = ic.b)
    end

    doubleoutput("Building simulation", label)
    simulation = Simulation(model, Δt = minimum([max_Δt/10, phys_params.T/100]), stop_time = duration, wall_time_limit = wall_time_limit)

    # ### The `TimeStepWizard`
    #
    # The TimeStepWizard manages the time-step adaptively, keeping the
    # Courant-Freidrichs-Lewy (CFL) number close to `1.0` while ensuring
    # the time-step does not increase beyond the maximum allowable value
    if :diffusive_cfl in keys(params)
        wizard = TimeStepWizard(cfl = 0.5, diffusive_cfl = params.diffusive_cfl, max_change = 1.1, max_Δt = max_Δt)
    else
        wizard = TimeStepWizard(cfl = 0.5, max_change = 1.1, max_Δt = max_Δt)
    end

    # Still some numerical noise at CFL 0.1 for Ri = 10⁴, but none for CFL = 0.05

    # A "Callback" pauses the simulation after a specified number of timesteps and calls a function (here the timestep wizard to update the timestep)
    # To update the timestep more or less often, change IterationInterval in the next line
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

    # ### A progress messenger
    # We add a callback that prints out a helpful progress message while the simulation runs.

    progress(sim) = begin
        string = @sprintf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e",
                            sim.model.clock.iteration,
                            prettytime(sim.model.clock.time),
                            prettytime(1e-9 * (time_ns() - start_time)),
                            sim.Δt,
                            AdvectiveCFL(sim.Δt)(sim.model))
        doubleoutput(string, label)
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

    # ### Output

    #=# Compute y-averages 𝐮̅(x,z) and b̅(x,z)
    u̅ = Field(Average(u, dims = 2))
    v̅ = Field(Average(v, dims = 2))
    w̅ = Field(Average(w, dims = 2))
    b̅ = Field(Average(b, dims = 2))
    ℬ = Field(w * b_pert)
    avg_ℬ = Field(Average(ℬ, dims = 2))=#

    doubleoutput("Setting Lagrangian particle output", label)
    filename = dir * "particles"
    simulation.output_writers[:particles] =
        JLD2Writer(model, (particles = model.particles,),
                                filename = filename * ".jld2",
                                schedule = TimeInterval(phys_params.T/30),
                                overwrite_existing = true)

    doubleoutput("Setting outptut for vertical slice", label)
    filename = dir * "BI_xz"
    simulation.output_writers[:xz_slices] =
        JLD2Writer(model, (; u, v, w, b, p, ζ, δ, fζ_g),
                                filename = filename * ".jld2",
                                indices = (:, 1, :),
                                schedule = TimeInterval(phys_params.T/30),
                                overwrite_existing = true)

    doubleoutput("Setting output for top slice", label)
    filename = dir * "BI_xy"
    simulation.output_writers[:xy_slices] =
        JLD2Writer(model, (; u, v, w, b, p, ζ, δ),
                                filename = filename * ".jld2",
                                indices = (:, :, resolution[3]),
                                schedule = TimeInterval(phys_params.T/30),
                                overwrite_existing = true,
                                with_halos = true)

    # Define a checkpointer, which will be activated when the simulation ends:
    simulation.output_writers[:checkpointer] = Checkpointer(model;
                dir = dir,
                schedule = WallTimeInterval(3wall_time_limit),
                prefix = "model_checkpoint",
                overwrite_existing = true)

    # Now, run the simulation
    doubleoutput("Simulation will last " * prettytime(duration), label)
    doubleoutput("Running simulation...\n", label)
    run!(simulation; checkpoint_at_end = true, pickup = pickup)

end

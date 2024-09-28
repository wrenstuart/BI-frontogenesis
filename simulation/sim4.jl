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

    u_x::Float64
    v_x::Float64
    u_z::Float64
    v_z::Float64
    w_x::Float64
    w_y::Float64
    b_x::Float64
    b_y::Float64
    b_z::Float64
    ζ::Float64
    δ::Float64

    fu_g::Float64
    fv_g::Float64
    fζ_g::Float64

    ∇ₕ²b_x::Float64
    ∇ₕ²b_y::Float64
    ∇ₕ²b_z::Float64
    b_xzz::Float64
    b_yzz::Float64
    b_zzz::Float64
    u_zxx::Float64
    u_zyy::Float64
    u_zzz::Float64
    v_zxx::Float64
    v_zyy::Float64
    v_zzz::Float64
    ζ_xx::Float64
    ζ_yy::Float64
    ζ_zz::Float64
    δ_xx::Float64
    δ_yy::Float64
    δ_zz::Float64
    ∇ₕ²ζ::Float64
    ∇ₕ²δ::Float64

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
    O() = zeros(n^2)
    if params.GPU
        x₀, y₀ = CuArray.([x₀, y₀])
        O() = CuArray(zeros(n^2))
    end
    particles = StructArray{MyParticle}((x₀, y₀, O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O(), O()))

    # Extract fundamental variable fields:
    velocities = VelocityFields(grid)
    u, v, w = velocities
    tracers = TracerFields((:b,), grid)
    b = tracers[1]
    pHY′ = CenterField(grid)
    pNHS = CenterField(grid)
    
    # Intermediary terms:
    p = pHY′ + pNHS
    u_y = ∂y(u)
    w_z = ∂z(w)

    # To store as auxillary fields:
    u_x = ∂x(u)
    v_x = ∂x(v)
    M²_on_f = phys_params.M²/f
    u_z = ∂z(u) - M²_on_f
    v_z = ∂z(v)
    w_x = ∂x(w)
    w_y = ∂y(w)
    b_x = ∂x(b)
    b_y = ∂y(b) + phys_params.M²
    b_z = ∂z(b)
    ζ = v_x - u_y
    δ = -w_z

    fu_g = -∂y(p)
    fv_g = ∂x(p)
    fζ_g = ∂x(fv_g) - ∂y(fu_g)  # Equivalently, ∇²Φ

    u_xxx = ∂x(∂x(u_x))
    u_xyy = ∂y(∂y(u_x))
    u_xzz = ∂z(∂z(u_x))
    v_xxx = ∂x(∂x(v_x))
    v_xyy = ∂y(∂y(v_x))
    v_xzz = ∂z(∂z(v_x))
    u_zxx = ∂x(∂x(u_z))
    u_zyy = ∂y(∂y(u_z))
    u_zzz = ∂z(∂z(u_z))
    v_zxx = ∂x(∂x(v_z))
    v_zyy = ∂y(∂y(v_z))
    v_zzz = ∂z(∂z(v_z))     # this one (specifically) is BADDDDDDD
    w_xxx = ∂x(∂x(w_x))
    w_xyy = ∂y(∂y(w_x))
    w_xzz = ∂z(∂z(w_x))
    w_yxx = ∂x(∂x(w_y))
    w_yyy = ∂y(∂y(w_y))
    w_yzz = ∂z(∂z(w_y))
    u_yxx = ∂x(∂x(u_y))
    u_yyy = ∂y(∂y(u_y))
    u_yzz = ∂z(∂z(u_y))
    ∇ₕ²(𝑓) = ∂x(∂x(𝑓)) + ∂y(∂y(𝑓))
    ζ_xx = ∂x(∂x(ζ))
    ζ_yy = ∂y(∂y(ζ))
    ζ_zz = ∂z(∂z(ζ))
    ∇ₕ²ζ = ∇ₕ²(ζ)
    δ_xx = ∂x(∂x(δ))
    δ_yy = ∂y(∂y(δ))
    δ_zz = ∂z(∂z(δ))
    ∇ₕ²δ = ∇ₕ²(δ)

    # The following do NOT work for tracked particles on a GPU:
    ∇ₕ²b_x = ∇ₕ²(b_x)
    b_xzz = ∂z(∂z(b_x))
    ∇ₕ²b_y = ∇ₕ²(b_y)
    b_yzz = ∂z(∂z(b_y))
    ∇ₕ²b_z = ∇ₕ²(b_z)
    b_zzz = ∂z(∂z(b_z))

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

    tracked_fields = (; ζ, δ, b_x, b_y, b_z, u_x, v_x, u_z, v_z, ∇ₕ²ζ, ζ_zz, ∇ₕ²δ, δ_zz, fζ_g, b_xzz, b_yzz)
    lagrangian_drifters = LagrangianParticles(particles; tracked_fields = tracked_fields)          

    # "Remember to use CuArray instead of regular Array when storing particle locations and properties on the GPU"?????

    # Build the model
    model = NonhydrostaticModel(; grid,
              advection = params.advection_scheme(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
              timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
              tracers = (; b),  # Set the name(s) of any tracers; here, b is buoyancy
              velocities = velocities,
              auxiliary_fields = (; ζ, δ, b_x, b_y, b_z, u_x, v_x, u_z, v_z, ∇ₕ²ζ, ζ_zz, ∇ₕ²δ, δ_zz, fζ_g, b_xzz, b_yzz),
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
    ζ = Field(model.auxiliary_fields.ζ)
    δ = Field(model.auxiliary_fields.δ)
    u_x = Field(model.auxiliary_fields.u_x)
    v_x = Field(model.auxiliary_fields.v_x)
    u_z = Field(model.auxiliary_fields.u_z)
    v_z = Field(model.auxiliary_fields.v_z)
    w_x = Field(∂x(w))
    w_y = Field(∂y(w))
    fu_g = Field(-∂y(p))
    fv_g = Field(∂x(p))
    fζ_g = Field(model.auxiliary_fields.fζ_g)
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
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, fu_g, fv_g, fζ_g),
                                filename = filename * ".jld2",
                                indices = (:, 1, :),
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    # Output the slice z = 0
    filename = "raw_data/" * label * "_BI_xy"
    simulation.output_writers[:xy_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, fu_g, fv_g, fζ_g),#, ∇ₕ²b_x, b_xzz, ∇ₕ²b_y, b_yzz),
                                filename = filename * ".jld2",
                                indices = (:, :, resolution[3]),
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    # Output the slice x = 0
    filename = "raw_data/" * label * "_BI_yz"
    simulation.output_writers[:yz_slices] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, fu_g, fv_g, fζ_g),
                                filename = filename * ".jld2",
                                indices = (1, :, :),
                                schedule = TimeInterval(phys_params.T/20),
                                overwrite_existing = true)

    # Output a horizontal slice in the middle (verticall speaking)
    filename = "raw_data/" * label * "_BI_xy_mid"
    simulation.output_writers[:xy_slices_mid] =
        JLD2OutputWriter(model, (; u, v, w, b, ζ, δ, u_x, v_x, u_z, v_z, w_x, w_y, b_x, b_y, b_z, fu_g, fv_g, fζ_g),
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
# # [Wind- and convection-driven mixing in an ocean surface boundary layer](@id gpu_example)
#
# This example simulates mixing by three-dimensional turbulence in an ocean surface
# boundary layer driven by atmospheric winds and convection. It demonstrates:
#
#   * How to set-up a grid with varying spacing in the vertical direction
#   * How to use the `SeawaterBuoyancy` model for buoyancy with a linear equation of state.
#   * How to use a turbulence closure for large eddy simulation.
#   * How to use a function to impose a boundary condition.
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add Oceananigans, JLD2, Plots"
# ```

# We start by importing all of the packages and functions that we'll need for this
# example.

using Random
using Printf
using Plots
using JLD2
using StructArrays

using Oceananigans
using Oceananigans.Units: minute, minutes, hour, day
using Oceananigans.Fields
using Oceananigans.Architectures: arch_array

#include("tracer_tendencies.jl")

# Define the grid
Nz = 24 # number of points in the vertical direction
Lz = 32 # domain depth
Nx = 64
Lx = 32
Ny = 64
Ly = 32

refinement = 1.2 # controls spacing near surface (higher means finer spaced)
stretching = 8   # controls rate of stretching at bottom 

## Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz

## Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement

## Bottom-intensified stretching function 
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

## Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

#grid = VerticallyStretchedRectilinearGrid(size = (128, 128, Nz), 
#                                          x = (0, 64),
#                                          y = (0, 64),
#                                          z_faces = z_faces)
#plot(grid.Δzᵃᵃᶜ[1:Nz], grid.zᵃᵃᶜ[1:Nz],
#     marker = :circle,
#     ylabel = "Depth (m)",
#     xlabel = "Vertical spacing (m)",
#     legend = nothing)

grid = RectilinearGrid(
                size=(Nx, Ny, Nz),
                x = (0, Lx),
                y = (0, Ly),
                z = (-Lz, 0))


# ## Buoyancy that depends on temperature and salinity
#
# We use the `SeawaterBuoyancy` model with a linear equation of state,

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion = 2e-4,
                                                                    haline_contraction = 8e-4))
# where $α$ and $β$ are the thermal expansion and haline contraction
# coefficients for temperature and salinity.
#
# ## Boundary conditions
#
# We calculate the surface temperature flux associated with surface heating of
# 200 W m⁻², reference density `ρₒ`, and heat capacity `cᴾ`,

Qʰ = 100  # W m⁻², surface _heat_ flux
ρₒ = 1026 # kg m⁻³, average density at the surface of the world ocean
cᴾ = 3991 # J K⁻¹ kg⁻¹, typical heat capacity for seawater

Qᵀ = Qʰ / (ρₒ * cᴾ) # K m s⁻¹, surface _temperature_ flux
nothing # hide

# Finally, we impose a temperature gradient `dTdz` both initially and at the
# bottom of the domain, culminating in the boundary conditions on temperature,

dTdz = 0 #0.01 # K m⁻¹


T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵀ),
                                bottom = GradientBoundaryCondition(dTdz))

# Note that a positive temperature flux at the surface of the ocean
# implies cooling. This is because a positive temperature flux implies
# that temperature is fluxed upwards, out of the ocean.
#
# For the velocity field, we imagine a wind blowing over the ocean surface
# with an average velocity at 10 meters `u₁₀`, and use a drag coefficient `cᴰ`
# to estimate the kinematic stress (that is, stress divided by density) exerted
# by the wind on the ocean:

u₁₀ = 0     # m s⁻¹, average wind velocity 10 meters above the ocean
cᴰ = 2.5e-3  # dimensionless drag coefficient
ρₐ = 1.225   # kg m⁻³, average density of air at sea-level

Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²

# The boundary conditions on `u` are thus

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))

# For salinity, `S`, we impose an evaporative flux of the form

@inline Qˢ(x, y, t, S, evaporation_rate) = - evaporation_rate * S # [salinity unit] m s⁻¹
nothing # hide

# where `S` is salinity. We use an evporation rate of 1 millimeter per hour,

evaporation_rate = 0e-3 / hour # m s⁻¹

# We build the `Flux` evaporation `BoundaryCondition` with the function `Qˢ`,
# indicating that `Qˢ` depends on salinity `S` and passing
# the parameter `evaporation_rate`,

evaporation_bc = FluxBoundaryCondition(Qˢ, field_dependencies=:S, parameters=evaporation_rate)

# The full salinity boundary conditions are

S_bcs = FieldBoundaryConditions(top=evaporation_bc)

# Configure Lagrangian particles
n_particles = 100

x₀ = Lx*rand(n_particles)
y₀ = Ly*rand(n_particles)
z₀ = -Lz/2*ones(n_particles)

struct TestParticle{T}
    x :: T
    y :: T
    z :: T
    T :: T
    S :: T
    C :: T
    Tmod :: T
end

#T = zeros(n_particles)
#S = zeros(n_particles)
#C = zeros(n_particles)
#Tmod = zeros(n_particles)

T = arch_array(CPU(), zeros(n_particles))
S = arch_array(CPU(), zeros(n_particles))
C = arch_array(CPU(), zeros(n_particles))
Tmod = arch_array(CPU(), zeros(n_particles))

particles = StructArray{TestParticle}((x₀,y₀,z₀,T,S,C,Tmod))

function particle_dynamics(lagrangian_particles,model,Δt)
    model.particles.properties.z.+=0.001*Δt
    model.particles.properties.z .= min.(model.particles.properties.z,0)  

    model.tracers.C[10,10,10] += 1
    
end  

tracers = TracerFields((:T, :S, :C),grid)
T, S, C = tracers

velocities = VelocityFields(grid)
u, v, w = velocities

model = NonhydrostaticModel(; grid, buoyancy,
                            advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            tracers = (:T, :S, :C),
                            velocities = velocities,
                            coriolis = FPlane(f=1e-4),
                            closure = ScalarDiffusivity(ν=1e-3, κ=1e-3),
                            boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs,))

tracked_fields = model.tracers
lagrangian_particles = LagrangianParticles(particles; tracked_fields, dynamics=particle_dynamics)

# Create a forcing function to include vertical 'slip' advection for buoyant or dense tracers

#tracer_slip_function(x, y, z, t, C, params) = -(params.wₛ * C)
#nothing # hide
#tracer_slip_parameters = (wₛ = 10/day)   # vertical 'slip' velocity in meters/day
# We tell `Forcing` that our 'tracer_slip' model depends
# on the tracer concentration `C` and the parameters,
#tracer_slip_dynamics = Forcing(tracer_slip_function, field_dependencies = :C,
#                            parameters = tracer_slip_parameters)

#dcdz = ComputedField(∂z(model.tracers.C))

growing_and_grazing(x, y, z, t, C, params) = -params.μ₀ * C * 0
nothing # hide

# with parameters

plankton_dynamics_parameters = (μ₀ = 1/day,   # surface growth rate
                                 λ = 5,       # sunlight attenuation length scale (m)
                                 m = 0.1/day) # mortality rate due to virus and zooplankton grazing

# We tell `Forcing` that our plankton model depends
# on the plankton concentration `P` and the chosen parameters,

plankton_dynamics = Forcing(growing_and_grazing, field_dependencies = :C,
                            parameters = plankton_dynamics_parameters)

# ## Model instantiation
#
# We fill in the final details of the model here: upwind-biased 5th-order
# advection for momentum and tracers, 3rd-order Runge-Kutta time-stepping,
# Coriolis forces, and the `AnisotropicMinimumDissipation` closure
# for large eddy simulation to model the effect of turbulent motions at
# scales smaller than the grid scale that we cannot explicitly resolve.

model = NonhydrostaticModel(; grid, buoyancy,
                            advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            tracers = (:T, :S, :C),
                            velocities = velocities,
                            coriolis = FPlane(f=1e-4),
                            closure = ScalarDiffusivity(ν=1e-3, κ=1e-3),
                            particles = lagrangian_particles,
                            forcing = (C=plankton_dynamics,),
                            boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs,))

# Notes:
#
# * To use the Smagorinsky-Lilly turbulence closure (with a constant model coefficient) rather than
#   `AnisotropicMinimumDissipation`, use `closure = SmagorinskyLilly()` in the model constructor.
#
# * To change the `architecture` to `GPU`, replace `architecture = CPU()` with
#   `architecture .= GPU()`.

# ## Initial conditions
#
# Our initial condition for temperature consists of a linear stratification superposed with
# random noise damped at the walls, while our initial condition for velocity consists
# only of random noise.

## Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

## Temperature initial condition: a stable density gradient with random noise superposed.
Tᵢ(x, y, z) = 20 + dTdz * z + 0 * model.grid.Lz * 1e-6 * Ξ(z)

## Velocity initial condition: random noise scaled by the friction velocity.
uᵢ(x, y, z) = 1e-3 * Ξ(z)
#uᵢ(x,y,z)=z/model.grid.Lz

Cᵢ(x,y,z)=exp(z/5)

## `set!` the `model` fields using functions or constants:
set!(model, u=uᵢ, w=uᵢ, T=Tᵢ, S=35, C=Cᵢ)

# ## Setting up a simulation
#
# We set-up a simulation with an initial time-step of 10 seconds
# that stops at 40 minutes, with adaptive time-stepping and progress printing.

simulation = Simulation(model, Δt=10.0, stop_time=10hour)

# The `TimeStepWizard` helps ensure stable time-stepping
# with a Courant-Freidrichs-Lewy (CFL) number of 1.0.

wizard = TimeStepWizard(cfl=1.0, max_change=1.1, max_Δt=1minute)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Nice progress messaging is helpful:

## Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.w),
                                prettytime(sim.run_wall_time))

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(10))

# We then set up the simulation:

# ## Output
#
# We use the `JLD2OutputWriter` to save ``x, z`` slices of the velocity fields,
# tracer fields, and eddy diffusivities. The `prefix` keyword argument
# to `JLD2OutputWriter` indicates that output will be saved in
# `ocean_wind_mixing_and_convection.jld2`.


# Diagnose the vertical vorticity
u, v, w = model.velocities # unpack velocity `Field`s
## Vertical vorticity [s⁻¹]
vorticity = ∂x(v) - ∂y(u)
# Create a NamedTuple for vorticity
vorticity_tuple = (; ω = vorticity)


## Create a NamedTuple with eddy viscosity
#eddy_viscosity = (; νₑ = model.diffusivity_fields.νₑ)


# Horizontal slice
#simulation.output_writers[:slices] =
#    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_viscosity, vorticity_tuple),
#                          prefix = "ocean_wind_mixing_and_convection",
#                     field_slicer = FieldSlicer(k=Int(grid.Nz)),
#                         schedule = TimeInterval(1minute),
#                            force = true)

# Vertical slice
simulation.output_writers[:xz_slices] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, vorticity_tuple),
                          filename = "slice_xz.jld2",
                          indices = (:, Int(grid.Ny/2), :),
                         schedule = TimeInterval(1minute),
                           overwrite_existing = true)

 simulation.output_writers[:xy_slices] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, vorticity_tuple),
                          filename = "slice_xy.jld2",
                         indices = (:, :, Int(grid.Nz)),
                         schedule = TimeInterval(1minute),
                            overwrite_existing = true)

simulation.output_writers[:particles] = JLD2OutputWriter(model, (particles=model.particles,), 
                            filename = "particles.jld2",
                          schedule = TimeInterval(1minute),
                             overwrite_existing = true)

# We're ready:

run!(simulation)

# ## Turbulence visualization
#
# We animate the data saved in `ocean_wind_mixing_and_convection.jld2`.
# We prepare for animating the flow by creating coordinate arrays,
# opening the file, building a vector of the iterations that we saved
# data at, and defining functions for computing colorbar limits:

## Coordinate arrays
xw, yw, zw = nodes(model.velocities.w)
xT, yT, zT = nodes(model.tracers.T)

## Open the file with our data
file_xz = jldopen(simulation.output_writers[:xz_slices].filepath)
file_xy = jldopen(simulation.output_writers[:xy_slices].filepath)
file_particles = jldopen(simulation.output_writers[:particles].filepath)

## Extract a vector of iterations
iterations = parse.(Int, keys(file_xy["timeseries/t"]))

""" Returns colorbar levels equispaced between `(-clim, clim)` and encompassing the extrema of `c`. """
function divergent_levels(c, clim, nlevels=21)
    cmax = maximum(abs, c)
    levels = clim > cmax ? range(-clim, stop=clim, length=nlevels) : range(-cmax, stop=cmax, length=nlevels)
    return (levels[1], levels[end]), levels
end

""" Returns colorbar levels equispaced between `clims` and encompassing the extrema of `c`."""
function sequential_levels(c, clims, nlevels=20)
    levels = range(clims[1], stop=clims[2], length=nlevels)
    cmin, cmax = minimum(c), maximum(c)
    cmin < clims[1] && (levels = vcat([cmin], levels))
    cmax > clims[2] && (levels = vcat(levels, [cmax]))
    return clims, levels
end
nothing # hide

# We start the animation at `t = 10minutes` since things are pretty boring till then:

times = [file_xy["timeseries/t/$iter"] for iter in iterations]

Tmod_save=zeros(n_particles,size(iterations[1:end])[1])
pT_save=zeros(n_particles,size(iterations[1:end])[1])

anim = @animate for (i, iter) in enumerate(iterations[1:end])

    @info "Drawing frame $i from iteration $iter..."

    t = file_xy["timeseries/t/$iter"]
    u = file_xy["timeseries/u/$iter"][:, :, 1]
    v = file_xy["timeseries/v/$iter"][:, :, 1]
    w = file_xy["timeseries/w/$iter"][:, :, 1]
    ω = file_xy["timeseries/ω/$iter"][:, :, 1]

    T = file_xz["timeseries/T/$iter"][:, 1, :]
    S = file_xz["timeseries/S/$iter"][:, 1, :]
    C = file_xz["timeseries/C/$iter"][:, 1, :]
    global particles = file_particles["timeseries/particles/$iter"]
    Tmod_save[:,i]=particles.Tmod
    pT_save[:,i]=particles.T

    ulims, ulevels = divergent_levels(u,0.05)
    wlims, wlevels = divergent_levels(w, 2e-2)
    ωlims,ωlevels = divergent_levels(ω,1e-2)
    Tlims, Tlevels = sequential_levels(T, (19.7, 19.99))
    Slims, Slevels = sequential_levels(S, (35, 35.005))
    Clims, Clevels = sequential_levels(C, (0, 1))

    kwargs_xz = (linewidth=0, xlabel="x (m)", ylabel="z (m)", aspectratio=1,
              xlims=(0, grid.Lx), ylims=(-grid.Lz, 0))
    kwargs_xy = (linewidth=0, xlabel="x (m)", ylabel="y (m)", aspectratio=1,
              xlims=(0, grid.Lx), ylims=(0, grid.Ly))

    w_plot = contourf(xw, yw, w'; color=:balance, clims=wlims, levels=wlevels, kwargs_xy...)
    scatter!(particles.x,particles.y, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    T_plot = contourf(xT, zT, T'; color=:thermal, clims=Tlims, levels=Tlevels, kwargs_xz...)
    scatter!(particles.x,particles.z, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    S_plot = contourf(xT, zT, S'; color=:haline,  clims=Slims, levels=Slevels, kwargs_xz...)
    scatter!(particles.x,particles.z, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    C_plot = contourf(xT, zT, C'; color=:haline, clims=Clims, levels=Clevels, kwargs_xz...)
    scatter!(particles.x,particles.z, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    ω_plot = contourf(xT,yT,ω'; color=:balance, clims=ωlims, level=ωlevels, kwargs_xy...)
    scatter!(particles.x,particles.y, legend=false, color=:black, markerstrokewidth=0, markersize=3)

    u_plot = contourf(xT, yT, u'; color=:balance, clims=ulims, levels=ulevels, kwargs_xy...)    
    scatter!(particles.x,particles.y, legend=false, color=:black, markerstrokewidth=0, markersize=3)
      

    w_title = @sprintf("vertical velocity (m s⁻¹), t = %s", prettytime(t))
    T_title = "temperature (ᵒC)"
    S_title = "salinity (g kg⁻¹)"
    ω_title = "vorticity (s⁻¹)"
    u_title = "x (m/s)"
    #ν_title = "eddy viscosity (m² s⁻¹)"

    ## Arrange the plots side-by-side.
    plot(w_plot, ω_plot, T_plot, C_plot, layout=(2, 2), size=(1200, 600),
         title=[w_title u_title T_title S_title])
    #plot(w_plot, size=(1200,1200), title=ω_title)

    #iter == iterations[end] && close(file)
end

mp4(anim, "ocean_wind_mixing_and_convection.mp4", fps = 20) # hide

close(file_xy)
close(file_xz)
close(file_particles)
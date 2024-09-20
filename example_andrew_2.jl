# This script simulates KH or  instability in 3D using Oceananigans and track particles

# Load some standard libraries that we will need
using Printf
using Oceananigans
using CUDA
using StructArrays
using Random

# First, we need to set some physical parameters for the simulation
# Set the domain size in non-dimensional coordinates

const Lx = 16.00  # size in the x-direction (16.00 form linstab)
const Lz = 15  # size in the vertical (z) direction
const Ly = 3  # size in the y-direction

# Set the grid size (see paper) (40min for 384*384*64)
const Nx = 100 # number of gridpoints in the x-direction
const Nz = 100 # number of gridpoints in the z-direction
const Ny = 20 # number of gridpoints in the y-direction

# Some timestepping parameters
const max_Δt = 0.05 # maximum allowable timestep 
const duration = 1000 # The non-dimensional duration of the simulation

# Set the Reynolds number (Re=Ul/ν)
const Re = 100
# Set the Prandtl number (Pr=ν/κ)
const Pr = 8

# Parameters for the initial condition:
const R=sqrt(8) # ratio of shear and density interface widths from Salehipour et al.   
const Ri₀=0.16*sqrt(8) # parameter from Salehipour et al
const S₀=1 # maximum shear
const N₀=sqrt(S₀^2*Ri₀) # maximum buoyancy frequency
const h=1 # shear layer width
const n_particles_1 = 101 # number of particles of equal buoyancy
const n_particles_2 = 101 # number of particles of equal buoyancy
const n_particles_3 = 101 # number of particles of equal buoyancy
const n_particles_4 = 101 # number of particles of equal buoyancy
const n_particles_5 = 101 # number of particles in the x direction
const n_particles_6 = 101 # number of particles in the x direction
const n_particles_7 = 101 # number of particles in the x direction
const n_particles_8 = 101 # number of particles in the x direction
const n_particles_9 = 101 # number of particles in the x direction


# Set the amplitude of the random perturbation (kick)
const kick = 0.005



# construct a rectilinear grid using an inbuilt Oceananigans function
# Here, we use periodic (cyclic) boundary conditions in x and y
grid = RectilinearGrid(GPU(), size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = (0, Lz), topology = (Periodic, Periodic, Bounded))

# No boundary conditions explicitly set - the BCs will default to free-slip and no-flux in z

# Now, define a 'model' where we specify the grid, advection scheme, bcs, and other settings
model = NonhydrostaticModel(; grid, 
              advection = WENO(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
            timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
                tracers = (:b),  # Set the name(s) of any tracers, here b is buoyancy
               buoyancy = Buoyancy(model=BuoyancyTracer()), # this tells the model that b will act as the buoyancy (and influence momentum) 
                closure = (ScalarDiffusivity(ν = 1 / Re, κ = 1 / (Re * Pr))),  # set a constant kinematic viscosity and diffusivty, here just 1/Re since we are solving the non-dimensional equations 
                coriolis = nothing # this line tells the mdoel not to include system rotation (no Coriolis acceleration)
)

u, v, w = model.velocities # unpack velocity `Field`s
b = model.tracers.b # extract the buoyancy


# Set initial conditions
# Here, we start with a tanh function for buoyancy and add a random perturbation to the velocity. 
# Set the random seed to get the same randomly perturbed velocity

Random.seed!(12341)
uᵢ(x, y, z) = S₀ * h * tanh( (z - Lz/2) / h) + kick * randn()

Random.seed!(56781)
vᵢ(x, y, z) = kick * randn()

Random.seed!(910111)
wᵢ(x, y, z) = kick * randn()

bᵢ(x, y, z) = N₀^2 * (h/R) * tanh( (z - Lz/2) / (h/R))

# Set the Lagrangian particles

struct LP{T}
  x :: T
  y :: T
  z :: T
  b :: T
end

# Set first set of particles

x₀_1 = CuArray((Lx / 2) * ones(n_particles_1)) 

y₀_1 = CuArray((Ly / 2) * ones(n_particles_1))

Initialb_1 = CuArray((-0.95 : 0.019 : 0.95).*(N₀^2 * (h/R)))

z₀_1 = (atanh.(Initialb_1 ./ (N₀^2 * (h/R)))) .* (h/R) .+ Lz/2

particles_1 = StructArray{LP}((x₀_1, y₀_1, z₀_1, Initialb_1))

# Set second set of particles

x₀_2 = CuArray((Lx * 0.75) * ones(n_particles_2)) 

y₀_2 = CuArray((Ly / 2) * ones(n_particles_2))

Initialb_2 = CuArray((-0.95 : 0.019 : 0.95).*(N₀^2 * (h/R)))

z₀_2 = (atanh.(Initialb_2 ./ (N₀^2 * (h/R)))) .* (h/R) .+ Lz/2

particles_2 = StructArray{LP}((x₀_2, y₀_2, z₀_2, Initialb_2))

# Set third set of particles

x₀_3 = CuArray(zeros(n_particles_3)) 

y₀_3 = CuArray((Ly / 2) * ones(n_particles_3))

Initialb_3 = CuArray((-0.95 : 0.019 : 0.95).*(N₀^2 * (h/R)))

z₀_3 = (atanh.(Initialb_3 ./ (N₀^2 * (h/R)))) .* (h/R) .+ Lz/2

particles_3 = StructArray{LP}((x₀_3, y₀_3, z₀_3, Initialb_3))

# Set fourth set of particles

x₀_4 = CuArray((Lx * 0.25) * ones(n_particles_4)) 

y₀_4 = CuArray((Ly / 2) * ones(n_particles_4))

Initialb_4 = CuArray((-0.95 : 0.019 : 0.95).*(N₀^2 * (h/R)))

z₀_4 = (atanh.(Initialb_4 ./ (N₀^2 * (h/R)))) .* (h/R) .+ Lz/2

particles_4 = StructArray{LP}((x₀_4, y₀_4, z₀_4, Initialb_4))

# Set fifth set of particles

x₀_5 = CuArray((0:0.01:1).*Lx) 

y₀_5 = CuArray((Ly / 2) * ones(n_particles_5))

z₀_5 = CuArray((Lz / 2) * ones(n_particles_5))

Initialb_5 = CuArray(zeros(n_particles_5))

particles_5 = StructArray{LP}((x₀_5, y₀_5, z₀_5, Initialb_5))

# Set sixth set of particles

x₀_6 = CuArray((0:0.01:1).*Lx) 

y₀_6 = CuArray((Ly / 2) * ones(n_particles_6))

z₀_6 = CuArray((Lz / 2 + 0.05) * ones(n_particles_6))

Initialb_6 = N₀^2 .* (h/R) .* tanh.( (z₀_6 .- Lz/2) ./ (h/R))

particles_6 = StructArray{LP}((x₀_6, y₀_6, z₀_6, Initialb_6))

# Set seventh set of particles

x₀_7 = CuArray((0:0.01:1).*Lx) 

y₀_7 = CuArray((Ly / 2) * ones(n_particles_7))

z₀_7 = CuArray((Lz / 2 + 0.2) * ones(n_particles_7))

Initialb_7 = N₀^2 .* (h/R) .* tanh.( (z₀_7 .- Lz/2) ./ (h/R))

particles_7 = StructArray{LP}((x₀_7, y₀_7, z₀_7, Initialb_7))

# Set eighth set of particles

x₀_8 = CuArray((0:0.01:1).*Lx) 

y₀_8 = CuArray((Ly / 2) * ones(n_particles_8))

z₀_8 = CuArray((Lz / 2 - 0.05) * ones(n_particles_8))

Initialb_8 = N₀^2 .* (h/R) .* tanh.( (z₀_8 .- Lz/2) ./ (h/R))

particles_8 = StructArray{LP}((x₀_8, y₀_8, z₀_8, Initialb_8))

# Set nineth set of particles

x₀_9 = CuArray((0:0.01:1).*Lx) 

y₀_9 = CuArray((Ly / 2) * ones(n_particles_9))

z₀_9 = CuArray((Lz / 2 - 0.2) * ones(n_particles_9))

Initialb_9 = N₀^2 .* (h/R) .* tanh.( (z₀_9 .- Lz/2) ./ (h/R))

particles_9 = StructArray{LP}((x₀_9, y₀_9, z₀_9, Initialb_9))

# Set total particle

particles_total = [particles_1; particles_2; particles_3; particles_4; particles_5; particles_6; particles_7; particles_8; particles_9]

lagrangian_particles = LagrangianParticles(particles_total, tracked_fields = model.tracers)

model = NonhydrostaticModel(; grid, 
              particles = lagrangian_particles,
              advection = WENO(),  # Specify the advection scheme.  Another good choice is WENO() which is more accurate but slower
            timestepper = :RungeKutta3, # Set the timestepping scheme, here 3rd order Runge-Kutta
                tracers = (; b),  # Set the name(s) of any tracers, here b is buoyancy
               buoyancy = Buoyancy(model=BuoyancyTracer()), # this tells the model that b will act as the buoyancy (and influence momentum) 
                closure = (ScalarDiffusivity(ν = 1 / Re, κ = 1 / (Re * Pr))),  # set a constant kinematic viscosity and diffusivty, here just 1/Re since we are solving the non-dimensional equations 
               coriolis = nothing # this line tells the mdoel not to include system rotation (no Coriolis acceleration)
)


# Send the initial conditions to the model to initialize the variables
set!(model, u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ)

# Now, we create a 'simulation' to run the model for a specified length of time
simulation = Simulation(model, Δt = max_Δt, stop_time = duration)

# ### The `TimeStepWizard`
wizard = TimeStepWizard(cfl = 0.85, max_change = 1.1, max_Δt = max_Δt)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
:w
# ### A progress messenger
start_time = time_ns()
progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        sim.model.clock.time,
                        prettytime(1e-9 * (time_ns() - start_time)),
                        sim.Δt,
                        AdvectiveCFL(sim.Δt)(sim.model))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

u, v, w = model.velocities # unpack velocity `Field`s
b = model.tracers.b # extract the buoyancy

#ω = Field(∂z((Center, Center, Center), u) - ∂x((Center, Center, Center), w)) # The spanwise vorticity
#compute!(ω)

#χ = Field((1 / (Re * Pr)) * (∂x((Center, Center, Center), b)^2 + ∂y((Center, Center, Center), b)^2 + ∂z((Center, Center, Center), b)^2)) # The dissipation rate of buoyancy variance
#compute!(χ)

#Inter1 = Field(∂x((Center, Center, Center), u)^2 + ∂y((Center, Center, Center), u)^2 + ∂z((Center, Center, Center), u)^2)
#Inter2 = Field(∂x((Center, Center, Center), v)^2 + ∂y((Center, Center, Center), v)^2 + ∂z((Center, Center, Center), v)^2)
#Inter3 = Field(∂x((Center, Center, Center), w)^2 + ∂y((Center, Center, Center), w)^2 + ∂z((Center, Center, Center), w)^2)
#ϵ = Field((1 / Re) * (Inter1 + Inter2 + Inter3)) # The dissipation rate of kinetic energy
#compute!(ϵ)

# ### Output
# Store slices

filename1 = "Primary_quantities_H"
simulation.output_writers[:xz_slices] =
    JLD2OutputWriter(model, (; b,), 
                             filename = filename1 * ".jld2",
                              indices = (:, trunc(Int,Ny/2), :),
                             schedule = TimeInterval(1),
                   overwrite_existing = true)


# Store particles

filename2 = "Primary_particles_H"
simulation.output_writers[:particles] = 
   JLD2OutputWriter(model, (particles = model.particles, ),
                             filename = filename2 * ".jld2", 
                             schedule = TimeInterval(1),
                   overwrite_existing = true)


# Store buoyancy dissipation 

#filename3 = "Buoyancy_dissipation"
#simulation.output_writers[:fields] =
#    JLD2OutputWriter(model, (; χ), 
#                             filename = filename3 * ".jld2",
#                              indices = (:, :, :),
#                             schedule = TimeInterval(1),
#                   overwrite_existing = true)



# Now, run the simulation
run!(simulation)


# After the simulation is finished, plot the results and save a movie

#include("plot_particle.jl")
#include("plot_overlay.jl")
#include("plot_statistics.jl")
#include("plot_pdf_ratio.jl")
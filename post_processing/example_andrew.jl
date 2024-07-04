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
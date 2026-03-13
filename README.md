To use, run main.jl

First, must create a small file, simulation/inputs/\[name\].jl, containing one function, sim_params. Here is an example:

function sim_params()
    Ri = 1                                  # Richardson number
    s = 1e4                                 # Stratification parameter, N²/f²
    ν_v = 1e-3                              # Vertical viscosity
    ν_h = 1e+1                              # Horizontal viscosity
    GPU = true                              # Set to false to run on CPU
    res = (512, 512, 64)                    # Number of gridpoints to be used
    advection_scheme = Centered
    horizontal_hyperviscosity = false       # If set to true, hyperviscosity will be used in the horizontal direction only
    short_duration = false                  # If set to true, simulation runs for 1/20th the time
    diffusive_cfl = 0.2
    return (; GPU, res, Ri, s, ν_h, ν_v, advection_scheme, horizontal_hyperviscosity, short_duration, diffusive_cfl)
end

To run the simulation, use julia to run main.jl with argument \[name\]:
    julia simulation/main.jl \[name\]

Raw output data will be stored in the raw_data/\[name\] folder

Package dependencies:

- CUDA
- CairoMakie
- FFTW
- JLD2
- Oceananigans v0.104 or higher (tested with v105.4)
- OffsetArrays
- StructArrays
- TickTock
- Unroll
- Dates
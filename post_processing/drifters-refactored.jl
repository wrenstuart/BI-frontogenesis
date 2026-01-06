using CairoMakie
using Oceananigans
using JLD2
using Printf
include("pp-io.jl")

f = 1e-4
ν_v = 1e-3
ν_h = 1e+1

per(i::Int64, N::Int64) :: Int64 = mod(i-1, N) + 1  # For periodic arrays

function rearrange_tracked_drifter_data(label::String)

    # Takes particle data from data[iter][var][drifter_num]
    # format to data[drifter_num][iter][var]
    in_name = data_dir(label) * "particles.jld2"
    out_name = data_dir(label) * "particles-rearranged.jld2"
    if ~isfile(in_name)
        throw("Could not find drifter data associated with label \"" * label * "\"")
    end
    in_file = jldopen(in_name)
    iterations = parse.(Int, keys(in_file["timeseries/t"]))
    t = [in_file["timeseries/t/$iter"] for iter in iterations]
    raw_keys = keys(in_file["timeseries/particles/0"])
    N = length(in_file["timeseries/particles/0"].x)
        # Number of particles
    ptcl_data = in_file["timeseries/particles"]
    drifters = [[(; (raw_keys .=> [ptcl_data["$iter"][var][n]
                    for var in raw_keys])...)
                    for iter in iterations]
                    for n = 1 : N]
    @save out_name t drifters

end

function extract_tracked_drifter_data(label)

    filename = data_dir(label) * "particles-rearranged.jld2"
    if ~isfile(filename)
        @info "Rearranging output data..."
        rearrange_tracked_drifter_data(label)
    end
    file = jldopen(filename)
    t = file["t"]
    drifters = file["drifters"]
    return t, drifters

end

function nearby_gridpoints(::Face, Δx::AbstractFloat, Nx::Int64, x₀::AbstractFloat) :: Tuple{Int64, Int64, Float64}

    # Finds gridpoint indices either side of x₀ on a periodic 1D Face grid
    i₋ = Int(floor(x₀/Δx)) + 1
    r = x₀/Δx - floor(x₀/Δx)
    i₊ = i₋ % Nx + 1
    return i₋, i₊, r
    
end

function nearby_gridpoints(::Center, Δx::AbstractFloat, Nx::Int64, x₀::AbstractFloat) :: Tuple{Int64, Int64, Float64}

    # Finds gridpoint indices either side of x₀ on a periodic 1D Center grid
    i₋ = Int(floor(x₀/Δx - 0.5)) + 1
    r = x₀/Δx - 0.5 - floor(x₀/Δx - 0.5)
    if i₋ == 0 i₋ = Nx end
    i₊ = i₋ % Nx + 1
    return i₋, i₊, r
    
end

function interpolate(var::String, (ℓx, ℓy), (x₀, y₀)::Tuple{Float64, Float64}, data::FileData, iter::Int64)

    # Interpolates a variable located on an (ℓx, ℓy) grid to a point (x₀, y₀)
    if ℓx == Center() x = data.xᶜ elseif ℓx == Face() x = data.xᶠ else throw("ℓx must be Center() or Face()") end
    if ℓy == Center() y = data.yᶜ elseif ℓy == Face() y = data.yᶠ else throw("ℓy must be Center() or Face()") end
    
    i₋, i₊, r₁ = nearby_gridpoints(ℓx, data.Δx, data.Nx, x₀)
    j₋, j₊, r₂ = nearby_gridpoints(ℓy, data.Δy, data.Ny, y₀)

    f₋₋ = data.file["timeseries/$var/$iter"][i₋, j₋, 1]
    f₋₊ = data.file["timeseries/$var/$iter"][i₋, j₊, 1]
    f₊₋ = data.file["timeseries/$var/$iter"][i₊, j₋, 1]
    f₊₊ = data.file["timeseries/$var/$iter"][i₊, j₊, 1]
    
    return (1 - r₁) * (1 - r₂) * f₋₋ + (1 - r₁) * r₂ * f₋₊ + r₁ * (1 - r₂) * f₊₋ + r₁ * r₂ * f₊₊

end

function extract_interpolated_drifter_data(eul_data::FileData, var::String, (ℓx, ℓy), x::Vector{<:AbstractFloat}, y::Vector{<:AbstractFloat}, drifter_t::Vector{<:AbstractFloat})
    
    iterations = eul_data.iterations
    eul_t = [eul_data.file["timeseries/t/$iter"] for iter in iterations]
    if drifter_t != eul_t throw("Eulerian and Lagrangian timeseries do not match") end
    var_interpolated = [interpolate(var, (ℓx, ℓy),
                                   (x[i], y[i]), eul_data, iter)
                      for (i, iter) in enumerate(iterations)]
    return var_interpolated

end

#=
TO DO:

Function for general extraction of tracked Lagrangian data (including t)
Function for general extraction of interpolated Lagrangian data (including t)

=#
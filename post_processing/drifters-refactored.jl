using CairoMakie
using Oceananigans
using JLD2
using Printf
include("pp-io.jl")

f = 1e-4
ν_v = 1e-3
ν_h = 1e+1

per(i::Int64, N::Int64) :: Int64 = mod(i-1, N) + 1  # For periodic arrays

label = "test_extra_visc_low_res"

struct FileData    # Useful struct for storing general data associated with a file

    file::JLD2.JLDFile{JLD2.MmapIO}
    xᶜ::Vector{Float64}
    xᶠ::Vector{Float64}
    yᶜ::Vector{Float64}
    yᶠ::Vector{Float64}
    zᶜ::Vector{Float64}
    zᶠ::Vector{Float64}
    Lx::Float64
    Ly::Float64
    Lz::Float64
    Nx::Int64
    Ny::Int64
    Nz::Int64
    Δx::Float64
    Δy::Float64
    Δz::Float64
    iterations::Vector{Int64}
    
end

function topdata(label::String) :: FileData

    filename_xy_top = data_dir(label) * "BI_xy.jld2"
    file = jldopen(filename_xy_top)
    xᶜ = file["grid/xᶜᵃᵃ"][4:end-3] # Equivalent to nodes(δ)[1] in other bits of code
    xᶠ = file["grid/xᶠᵃᵃ"][4:end-3] # Equivalent to nodes(ζ)[1] in other bits of code
    yᶜ = file["grid/yᵃᶜᵃ"][4:end-3]
    yᶠ = file["grid/yᵃᶠᵃ"][4:end-3]
    zᶜ = file["grid/zᵃᵃᶜ"][4:end-3]
    zᶠ = file["grid/zᵃᵃᶠ"][4:end-3]
    Nx = length(xᶜ) # (xᶠ would equally work)
    Ny = length(yᶜ)
    Nz = length(zᶜ)
    Δx, Δy, Δz = (x -> x[2] - x[1]).([xᶜ, yᶜ, zᶜ])
    Lx = Δx * Nx
    Ly = Δy * Ny
    Lz = Δz * Nz
    iterations = parse.(Int64, keys(file["timeseries/t"]))
    return FileData(file, xᶜ, xᶠ, yᶜ, yᶠ, zᶜ, zᶠ, Lx, Ly, Lz, Nx, Ny, Nz, Δx, Δy, Δz, iterations)

end

function nearby_gridpoints(::Face, Δx::Float64, Nx::Int64, x₀::Float64) :: Tuple{Int64, Int64, Float64}

    # Finds gridpoint indices either side of x₀ on a periodic 1D Face grid
    i₋ = Int(floor(x₀/Δx)) + 1
    r = x₀/Δx - floor(x₀/Δx)
    i₊ = i₋ % Nx + 1
    return i₋, i₊, r
    
end

function nearby_gridpoints(::Center, Δx::Float64, Nx::Int64, x₀::Float64) :: Tuple{Int64, Int64, Float64}

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

function ζ_tseriess(label::String)
    
    # WE ASSUME THAT BI_XY AND PARTICLE ITERATIONS ARE THE SAME

    check_pp_dir(label)
    eul_data = topdata(label)
    ptcl_file = jldopen(data_dir(label) * "particles.jld2")
    iterations = eul_data.iterations

    t = [ptcl_file["timeseries/t/$iter"] for iter in iterations]
    x_tracked = [ptcl_file["timeseries/particles/$iter"].x[1] for iter in iterations]
    y_tracked = [ptcl_file["timeseries/particles/$iter"].y[1] for iter in iterations]
    ζ_tracked = [ptcl_file["timeseries/particles/$iter"].ζ[1] for iter in iterations]
    ζ_interpolated = [interpolate("ζ", (Face(), Face()),
                                   (x_tracked[i], y_tracked[i]), eul_data, iter)
                      for (i, iter) in enumerate(iterations)]

    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(f*t, ζ_tracked)
    lines!(f*t, ζ_interpolated)
    display(fig)

end
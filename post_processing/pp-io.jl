data_dir(label::String) :: String = "raw_data/" * label * "/"
pp_dir(label::String) :: String = "pretty_things/" * label * "/"
using JLD2

function check_pp_dir(label::String)
    if !isdir(pp_dir(label))
        mkdir(pp_dir(label))
    end
end

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
    zᶜ = file["grid/z/cᵃᵃᶜ"][4:end-3]
    zᶠ = file["grid/z/cᵃᵃᶠ"][4:end-3]
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
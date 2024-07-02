#=using CairoMakie
using Oceananigans
using JLD2
using Printf
using Oceananigans.Units
using Makie.Colors
using OffsetArrays=#

using CairoMakie
using Oceananigans
using JLD2
using Printf

struct FileData
    file::JLD2.JLDFile{JLD2.MmapIO}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    Lx::Float64
    Ly::Float64
    Lz::Float64
    Nx::Int64
    Ny::Int64
    Nz::Int64
    Δx::Float64
    Δy::Float64
    Δz::Float64
end

function FileData(file::JLD2.JLDFile{JLD2.MmapIO}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64})
    Nx, Ny, Nz = length.([x, y, z])
    Lx = x[Nx] - x[1]
    Ly = y[Ny] - y[1]
    Lz = z[Nz] - z[1]
    FileData(file, x, y, z, Lx, Ly, Lz, Nx, Ny, Nz, Lx/Nx, Ly/Ny, Lz/Nz)
end

function topdata(label::String)
    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"
    file = jldopen(filename_xy_top)
    grid_x, grid_y, grid_z = nodes(FieldTimeSeries(filename_xy_top, "ζ₃", iterations = 0))
    data = FileData(file, [x for x in grid_x], [y for y in grid_y], [z for z in grid_z])
    return data
end

function grid_interpolate(data::FileData, var::String, x::Float64, y::Float64, iter::Int)

    i₋ = Int(floor(x/data.Lx * data.Nx)) + 1
    i₊ = i₋ % data.Nx + 1
    j₋ = Int(floor(y/data.Ly * data.Ny)) + 1
    j₊ = j₋ % data.Ny + 1
    x_frac = (x - data.x[i₋]) / data.Δx
    y_frac = (y - data.y[j₋]) / data.Δy

    f₋₋ = data.file["timeseries/$var/$iter"][i₋, j₋, 1]
    f₋₊ = data.file["timeseries/$var/$iter"][i₋, j₊, 1]
    f₊₋ = data.file["timeseries/$var/$iter"][i₊, j₋, 1]
    f₊₊ = data.file["timeseries/$var/$iter"][i₊, j₊, 1]
    f = x_frac * y_frac * f₋₋ + x_frac * (1-y_frac) * f₋₊ + (1-x_frac) * y_frac * f₊₋ + (1-x_frac) * (1-y_frac) * f₊₊
    
    return f

end

function tracer_release(label::String, 𝐱₀::Vector{Float64})

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"
    file = jldopen(filename_xy_top)
    grid_x, grid_y, grid_z = nodes(FieldTimeSeries(filename_xy_top, "ζ₃", iterations = 0))
    data = FileData(file, [x for x in grid_x], [y for y in grid_y], [z for z in grid_z])

    iterations = parse.(Int, keys(file["timeseries/t"]))
    t = [file["timeseries/t/$iter"] for iter in iterations]

    𝐱 = [[0.0, 0.0] for i = 1: length(iterations)]
    𝐱[1] = 𝐱₀

    for n = 1 : length(iterations) - 1

        iter = iterations[n]
        𝐱₁ = 𝐱[n]
        𝐮₁ = [grid_interpolate(data, "u", 𝐱₁[1], 𝐱₁[2], iter), grid_interpolate(data, "v", 𝐱₁[1], 𝐱₁[2], iter)]
        𝐱₂ = (𝐱[n] + (t[n+1]-t[n]) * 𝐮₁ + [data.Lx, data.Ly]) .% [data.Lx, data.Ly]
        𝐮₂ = [grid_interpolate(data, "u", 𝐱₂[1], 𝐱₂[2], iter), grid_interpolate(data, "v", 𝐱₂[1], 𝐱₂[2], iter)]
        𝐮 = (𝐮₁ + 𝐮₂) / 2
        𝐱[n+1] = (𝐱[n] + (t[n+1]-t[n]) * 𝐮 + [data.Lx, data.Ly]) .% [data.Lx, data.Ly]

    end

    return 𝐱

end

function tracer_release(data::FileData, 𝐱₀::Vector{Float64})

    iterations = parse.(Int, keys(data.file["timeseries/t"]))
    t = [data.file["timeseries/t/$iter"] for iter in iterations]

    𝐱 = [[0.0, 0.0] for i = 1: length(iterations)]
    𝐱[1] = 𝐱₀

    for n = 1 : length(iterations) - 1

        iter = iterations[n]
        𝐱₁ = 𝐱[n]
        𝐮₁ = [grid_interpolate(data, "u", 𝐱₁[1], 𝐱₁[2], iter), grid_interpolate(data, "v", 𝐱₁[1], 𝐱₁[2], iter)]
        𝐱₂ = (𝐱[n] + (t[n+1]-t[n]) * 𝐮₁ + [data.Lx, data.Ly]) .% [data.Lx, data.Ly]
        𝐮₂ = [grid_interpolate(data, "u", 𝐱₂[1], 𝐱₂[2], iter), grid_interpolate(data, "v", 𝐱₂[1], 𝐱₂[2], iter)]
        𝐮 = (𝐮₁ + 𝐮₂) / 2
        𝐱[n+1] = (𝐱[n] + (t[n+1]-t[n]) * 𝐮 + [data.Lx, data.Ly]) .% [data.Lx, data.Ly]

    end

    return 𝐱

end

function tracer_grid(data::FileData, n::Int)
    # create an equally spaced n×n grid of tracers in the initial conditions
    return vec([tracer_release(data, [data.Lx*(i-1)/n, data.Ly*(j-1)/n]) for i = 1:n, j = 1:n])
end

function extract_tracers(label::String)
    
    filename_tracers = "raw_data/" * label * "_particles.jld2"
    file = jldopen(filename_tracers)
    tracers_prim = [[file["timeseries/particles/$iter"].x, file["timeseries/particles/$iter"].y] for iter in parse.(Int, keys(file["timeseries/t"]))]
    # Above is indexed by [iter][x/y][tracer_number]
    n_iters = length(tracers_prim)
    n_tracers = length(tracers_prim[1][1])
    return [[[tracers_prim[i][j][k] for j = 1 : 2] for i = 1 : n_iters] for k = 1 : n_tracers]

end

function lagr_track(data::FileData, var_label::String, drifter::Vector{Vector{Float64}})
    return lagr_func_track(data, x -> x[1], [var_label], drifter)
end

function lagr_func_track(data::FileData, f::Function, input_labels::Vector{String}, drifter::Vector{Vector{Float64}})

    # f is a (scalar) function on a tuple X, the value of which are determined by input_labels
    # E.g. one could have f = 𝐮 -> (𝐮[1]^2 + 𝐮[2]^2 + 𝐮[3]^2) / 2 for kinetic energy, and
    # input_labels would be ["u", "v", "w"] here
    # (note that each input of f must be treated as a vector)

    iterations = parse.(Int, keys(data.file["timeseries/t"]))
    t = [data.file["timeseries/t/$iter"] for iter in iterations]
    output = zeros(Float64, length(iterations))
    for i = 1 : length(iterations)
        iter = iterations[i]
        t[i] = data.file["timeseries/t/$iter"]
        x, y = drifter[i]
        input = [grid_interpolate(data, var, x, y, iter) for var in input_labels]
        output[i] = f(input)
    end

    return t, output

end

function tracer_δ(label)

    f = 1e-4
    data = topdata(label)
    drifters = tracer_grid(data, 3)[1:5]
    fig = Figure()
    ax = Axis(fig[1, 1])
    for drifter in drifters
        t, δ = lagr_func_track(data, x -> x[1]/f, ["δ"], drifter)
        i₁ = Int(round(length(t)/3))
        i₂ = Int(round(2length(t)/3))
        δ_smooth = zeros(length(δ))
        for i = i₁ : i₂
            δ_smooth[i] = sum(δ[i-4:i+4])/9
        end
        lines!(ax, t[i₁:i₂], δ_smooth[i₁:i₂])
    end
    display(fig)
    save("pretty_things/tracer-delta_" * label * ".pdf", fig)
    
end

function ani_tracers(label::String)

    data = topdata(label)
    iterations = parse.(Int, keys(data.file["timeseries/t"]))
    tracers = tracer_grid(data, 5)
    f = 1e-4

    frame = Observable(1)

    tracers_now_x = lift(i -> [tracer[i][1]/1e3 for tracer in tracers], frame)
    tracers_now_y = lift(i -> [tracer[i][2]/1e3 for tracer in tracers], frame)

    ζ_on_f = lift(frame) do i
        iter = iterations[i]
        data.file["timeseries/ζ₃/$iter"][:, :, 1]/f
    end

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1)
    heatmap!(ax, data.x/1e3, data.y/1e3, ζ_on_f, colormap = :coolwarm, colorrange = (-20, 20));
    scatter!(ax, tracers_now_x, tracers_now_y, marker = '.', markersize = 30, color = :black)

    record(i -> frame[] = i, fig, "pretty_things/tracer_" * label * ".mp4", 1 : length(iterations), framerate = 20)

end

function tracer_δ_2(label)

    f = 1e-4
    data = topdata(label)
    drifters = extract_tracers(label)[1:5]
    fig = Figure()
    ax = Axis(fig[1, 1])
    for drifter in drifters
        t, δ = lagr_func_track(data, x -> x[1]/f, ["δ"], drifter)
        i₁ = Int(round(length(t)/3))
        i₂ = Int(round(2length(t)/3))
        δ_smooth = zeros(length(δ))
        for i = i₁ : i₂
            δ_smooth[i] = sum(δ[i-4:i+4])/9
        end
        lines!(ax, t[i₁:i₂], δ_smooth[i₁:i₂])
    end
    display(fig)
    save("pretty_things/tracer-delta_" * label * ".pdf", fig)
    
end

function ani_tracers_2(label::String)
    
    data = topdata(label)
    iterations = parse.(Int, keys(data.file["timeseries/t"]))
    tracers = extract_tracers(label)
    f = 1e-4

    frame = Observable(1)

    tracers_now_x = lift(i -> [tracer[i][1]/1e3 for tracer in tracers], frame)
    tracers_now_y = lift(i -> [tracer[i][2]/1e3 for tracer in tracers], frame)

    ζ_on_f = lift(frame) do i
        iter = iterations[i]
        data.file["timeseries/ζ₃/$iter"][:, :, 1]/f
    end

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1)
    heatmap!(ax, data.x/1e3, data.y/1e3, ζ_on_f, colormap = :coolwarm, colorrange = (-20, 20));
    scatter!(ax, tracers_now_x, tracers_now_y, marker = '.', markersize = 30, color = :black)

    record(i -> frame[] = i, fig, "pretty_things/tracer_" * label * ".mp4", 1 : length(iterations), framerate = 20)
    
end
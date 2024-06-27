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

function grid_interpolate(file::JLD2.JLDFile{JLD2.MmapIO}, grid::NamedTuple, var::String, x::Float64, y::Float64, iter::Int)

    i₋ = Int(floor(x/grid.Lx * grid.N)) + 1
    i₊ = i₋ % grid.N + 1
    j₋ = Int(floor(y/grid.Ly * grid.M)) + 1
    j₊ = j₋ % grid.M + 1
    x_frac = (x - grid.x[i₋]) / grid.Δx
    y_frac = (y - grid.y[j₋]) / grid.Δy

    f₋₋ = file["timeseries/$var/$iter"][i₋, j₋, 1]
    f₋₊ = file["timeseries/$var/$iter"][i₋, j₊, 1]
    f₊₋ = file["timeseries/$var/$iter"][i₊, j₋, 1]
    f₊₊ = file["timeseries/$var/$iter"][i₊, j₊, 1]
    f = x_frac * y_frac * f₋₋ + x_frac * (1-y_frac) * f₋₊ + (1-x_frac) * y_frac * f₊₋ + (1-x_frac) * (1-y_frac) * f₊₊
    
    return f

end

function tracer_release(label::String, 𝐱₀::Vector{Float64})

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"
    file = jldopen(filename_xy_top)
    iterations = parse.(Int, keys(file["timeseries/t"]))
    t = [file["timeseries/t/$iter"] for iter in iterations]

    u_ic = FieldTimeSeries(filename_xy_top, "ζ₃", iterations = 0)
    grid_x, grid_y, ~ = nodes(u_ic)
    Lx = grid_x[end]
    Ly = grid_y[end]
    N = length(grid_x)
    M = length(grid_y)
    Δx = Lx / N
    Δy = Ly / M
    grid = (x = grid_x, y = grid_y, Lx = Lx, Ly = Ly, M = M, N = N, Δx = Δx, Δy = Δy)

    𝐱 = [[0.0, 0.0] for i = 1: length(iterations)]
    𝐱[1] = 𝐱₀

    for n = 1 : length(iterations) - 1

        iter = iterations[n]
        𝐱₁ = 𝐱[n]
        𝐮₁ = [grid_interpolate(file, grid, "u", 𝐱₁[1], 𝐱₁[2], iter), grid_interpolate(file, grid, "v", 𝐱₁[1], 𝐱₁[2], iter)]
        𝐱₂ = (𝐱[n] + (t[n+1]-t[n]) * 𝐮₁ + [Lx, Ly]) .% [Lx, Ly]
        𝐮₂ = [grid_interpolate(file, grid, "u", 𝐱₂[1], 𝐱₂[2], iter), grid_interpolate(file, grid, "v", 𝐱₂[1], 𝐱₂[2], iter)]
        𝐮 = (𝐮₁ + 𝐮₂) / 2
        𝐱[n+1] = (𝐱[n] + (t[n+1]-t[n]) * 𝐮 + [Lx, Ly]) .% [Lx, Ly]

    end

    return 𝐱

end

function tracer_grid(label::String, n::Int)

    # create an equally spaced n×n grid of tracers in the initial conditions

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"
    u_ic = FieldTimeSeries(filename_xy_top, "ζ₃", iterations = 0)
    grid_x, grid_y, ~ = nodes(u_ic)
    Lx = grid_x[end]
    Ly = grid_y[end]

    return vec([tracer_release(label, [Lx*(i-1)/n, Ly*(j-1)/n]) for i = 1:n, j = 1:n])

end

function ani_tracers(label::String)

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"
    ζ_ic = FieldTimeSeries(filename_xy_top, "ζ₃", iterations = 0)
    x_grid, y_grid, z_grid = nodes(ζ_ic)
    file = jldopen(filename_xy_top)
    iterations = parse.(Int, keys(file["timeseries/t"]))

    tracers = tracer_grid(label, 5)
    f = 1e-4

    frame = Observable(1)

    tracers_now_x = lift(i -> [tracer[i][1]/1e3 for tracer in tracers], frame)
    tracers_now_y = lift(i -> [tracer[i][2]/1e3 for tracer in tracers], frame)

    ζ_on_f = lift(frame) do i
        iter = iterations[i]
        file["timeseries/ζ₃/$iter"][:, :, 1]/f
    end

    fig = Figure()
    ax = Axis(fig[1, 1])
    heatmap!(ax, x_grid/1e3, y_grid/1e3, ζ_on_f, colormap = :coolwarm, colorrange = (-20, 20));
    scatter!(ax, tracers_now_x, tracers_now_y, marker = '.', markersize = 30, color = :black)

    record(i -> frame[] = i, fig, "pretty_things/tracer_" * label * ".mp4", Int64(round(length(iterations)*0.5)) : length(iterations), framerate = 20)

end
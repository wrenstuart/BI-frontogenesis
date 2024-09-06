# Old functions for releasing drifters in post_processing (removed in favour of
# drifters handled by OCeananigans during simulation)

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
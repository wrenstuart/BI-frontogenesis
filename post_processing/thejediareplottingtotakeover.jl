using Oceananigans, JLD2, Makie, Printf
using Oceananigans.Units
using CairoMakie



function surface_stats(file_label, var_label)
    f = x -> x[1]
    input_labels = [var_label]
    return surface_function_stats(file_label, f, input_labels)
end

function surface_function_stats(file_label, f, input_labels)
    # f is a function on a tuple X, the value of which are determined by input_labels
    # E.g. one could have f = 𝐮 -> (𝐮[1].^2 + 𝐮[2].^2 + 𝐮[3].^2) / 2 for kinetic energy,
    # and input_labels would be ["u", "v", "w"] here
    # (note that each input of f must be treated as a vector)
    filename = "raw_data/" * file_label * "_BI_xy.jld2"
    @info filename
    @info pwd()
    file = jldopen(filename)
    ic = FieldTimeSeries(filename, input_labels[1], iterations = 0)
    x, y, z = nodes(ic)
    iterations = parse.(Int, keys(file["timeseries/t"]))
    t = zeros(length(iterations))   # times that the iterations represent
    μ = zeros(Float64, length(iterations))
    σ² = zeros(Float64, length(iterations))
    skew = zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t[i] = file["timeseries/t/$iter"]
        X = [file["timeseries/" * var_label * "/$iter"] for var_label in input_labels]
        N = Float64(length(x) * length(y))
        fX = f(X)
        μ[i] = sum(fX) / N
        σ²[i] = sum((fX .- μ[i]) .^ 2) / N
        skew[i] = sum(((fX .- μ[i]) .^ 3) / (σ²[i] .^ 1.5)) / N
    end

    return t[11:end], μ[11:end], σ²[11:end], skew[11:end]
end



#f = 𝐮 -> (𝐮[1].^2 + 𝐮[2].^2 + 𝐮[3].^2) / 2
#input_labels = ["u", "v", "w"]
#t, μ, σ², skew = surface_function_stats("test", f, input_labels)

t, μ_ζ, σ²_ζ, skew_ζ = surface_stats("test", "ζ₃")
t, μ_δ, σ²_δ, skew_δ = surface_stats("test", "δ")

lines(t .* 1e-4, skew_δ)
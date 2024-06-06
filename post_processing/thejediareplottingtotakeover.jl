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

function buoyancy_flux(file_label)
    filename = "raw_data/" * file_label * "_BI_y-avg.jld2"
    @info filename
    @info pwd()
    file = jldopen(filename)
    ic = FieldTimeSeries(filename, "avg_ℬ", iterations = 0)
    x, y, z = nodes(ic)
    iterations = parse.(Int, keys(file["timeseries/t"]))
    t = zeros(length(iterations))   # times that the iterations represent
    μ = zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t[i] = file["timeseries/t/$iter"]
        X = file["timeseries/avg_ℬ/$iter"]
        N = Float64(length(x) * length(z))
        μ[i] = sum(X) / N
    end

    return t[11:end], μ[11:end]
end

label = "test"
function basic_plots(label, Ri)
    f = 1e-4
    s = 1e4
    H = 50
    U = (s/Ri)^0.5 * H * f
    t, 𝒦, ~, ~ = surface_function_stats(label, 𝐮 -> ((𝐮[1] .- U).^2 + 𝐮[2].^2 + 𝐮[3].^2) / 2, ["u", "v", "w"])
    t, K, ~, ~ = surface_function_stats(label, 𝐮 -> (𝐮[1].^2 + 𝐮[2].^2 + 𝐮[3].^2) / 2, ["u", "v", "w"])
    t, ℬ_surf, ~, ~ = surface_function_stats(label, x -> - x[1] .* x[2], ["w", "b"])
    #t, ℬ_vol = buoyancy_flux(label)
    ℬ = ℬ_surf
    t, μ_ζ, σ²_ζ, skew_ζ = surface_stats(label, "ζ₃")
    t, μ_δ, σ²_δ, skew_δ = surface_stats(label, "δ")
    #t, μ_Ri, ~, ~ = surface_function_stats("test", x -> - f^2 .* x[1] ./ (u_z.^2 .+ v_z.^2), ["b_z", "u_z", "v_z"])
    μ_Ri = t .* 0

    fig, ax, l1 = lines(t .* 1e-4, skew_ζ)
    l2 = lines!(t .* 1e-4, skew_δ)
    save("pretty_things/" * label * "_skew.pdf", fig)
    fig, ax, l1 = lines(t .* 1e-4, log.(σ²_ζ.^0.5/f))
    l2 = lines!(t .* 1e-4, log.(σ²_δ.^0.5/f))
    save("pretty_things/" * label * "_var.pdf", fig)
    fig, ax, l = lines(t .* 1e-4, μ_Ri)
    save("pretty_things/" * label * "_Ri.pdf", fig)
    fig, ax, l = lines(t .* 1e-4, ℬ)
    save("pretty_things/" * label * "_ℬ.pdf", fig)
end
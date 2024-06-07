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
function timeseries(label, Ri)

    f = 1e-4
    s = 1e4
    H = 50
    U = (s/Ri)^0.5 * H * f
    t, 𝒦, ~, ~ = surface_function_stats(label, 𝐮 -> ((𝐮[1] .- U).^2 + 𝐮[2].^2 + 𝐮[3].^2) / 2, ["u", "v", "w"])
    t, K, ~, ~ = surface_function_stats(label, 𝐮 -> (𝐮[1].^2 + 𝐮[2].^2 + 𝐮[3].^2) / 2, ["u", "v", "w"])
    t, ℬ_surf, ~, ~ = surface_function_stats(label, x -> - x[1] .* x[2], ["w", "b"])
    t, ℬ_vol = buoyancy_flux(label)
    ℬ = ℬ_vol
    t, μ_ζ, σ²_ζ, skew_ζ = surface_stats(label, "ζ₃")
    t, μ_δ, σ²_δ, skew_δ = surface_stats(label, "δ")
    #t, μ_Ri, ~, ~ = surface_function_stats(label, x -> - f^2 .* x[1] ./ (x[2].^2 .+ x[3].^2), ["b_z", "u_z", "v_z"])
    μ_Ri = t .* 0

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"ft", ylabel=L"\text{Skewness}")
    l1 = lines!(t .* 1e-4, skew_ζ, label = L"\zeta")
    l2 = lines!(t .* 1e-4, skew_δ, label = L"\delta")
    axislegend(position = :rc)
    save("pretty_things/" * label * "_skew.pdf", fig)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"ft", ylabel=L"\text{log-RMS}")
    l1 = lines!(t .* 1e-4, log.(σ²_ζ.^0.5/f), label = L"\zeta")
    l2 = lines!(t .* 1e-4, log.(σ²_δ.^0.5/f), label = L"\delta")
    axislegend(position = :rc)
    save("pretty_things/" * label * "_var.pdf", fig)
    #fig, ax, l = lines(t .* 1e-4, μ_Ri)
    #save("pretty_things/" * label * "_Ri.pdf", fig)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"ft", ylabel=L"\text{Buoyancy flux }\mathcal{B}=-\langle wb\rangle")
    l = lines!(t .* 1e-4, ℬ)
    save("pretty_things/" * label * "_ℬ.pdf", fig)

end

function get_iter(label, t)
    filename = "raw_data/" * label * "_BI_xy.jld2"
    file = jldopen(filename)
    iterations = parse.(Int, keys(file["timeseries/t"]))
    ts = [file["timeseries/t/$iter"] for iter in iterations]
    i = argmin(abs.(ts .- t))
    return iterations[i]
end

function ζ_δ_joint_freq(ζ, δ)

    f = 1e-4
    ζ = vec(ζ)
    δ = vec(δ)
    n = 200
    a = -12f
    b = 12f
    Δ  = (b-a) / (n-2)
    midpoints = a + 0.5Δ : Δ : b - 0.5Δ
    freq = zeros(Float64, (n, n))
    
    for k = 1:length(ζ)
        if ζ[k] < b && ζ[k] > a && δ[k] < b && δ[k] > a
            i = argmin(abs.(ζ[k] .- midpoints))
            j = argmin(abs.(δ[k] .- midpoints))
            freq[i, j] += 1
        end
    end

    return midpoints, midpoints, freq

end

function snapshots(label, t)

    f = 1e-4
    iter = get_iter(label, t)
    filename = "raw_data/" * label * "_BI_xy.jld2"
    file = jldopen(filename)
    ic = FieldTimeSeries(filename, "ζ₃", iterations = 0)
    x, y, z = nodes(ic)
    ζ = file["timeseries/ζ₃/$iter"][:, :, 1]
    δ = file["timeseries/δ/$iter"][:, :, 1]
    b = file["timeseries/b/$iter"][:, :, 1]

    fig = Figure()
    ax = Axis(fig[1, 1], ylabel=L"\text{Probability density}")
    h1 = stephist!(vec(ζ/f), bins = 1000, normalization = :pdf, label = L"\zeta/f")
    h2 = stephist!(vec(δ/f), bins = 1000, normalization = :pdf, label = L"\delta/f")
    axislegend()
    xlims!(ax, (-6, 6))
    display(fig)
    save("pretty_things/" * label * "_hist.pdf", fig)

    ζ_ax, δ_ax, freq = ζ_δ_joint_freq(ζ, δ)
    fig = Figure(resolution = (600,600), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = L"\zeta/f", ylabel = L"\delta/f",
    title = L"\text{Joint histogram of vorticity and divergence}")
    hm = heatmap!(ax, ζ_ax/f, δ_ax/f, -log.(1 .+ freq), colormap = :bilbao)
    lines!([0, 0], [δ_ax[1]/f, δ_ax[end]/f], color = :black, linestyle = :dash)
    lines!([ζ_ax[1]/f, ζ_ax[end]/f], [0, 0], color = :black, linestyle = :dash)
    save("pretty_things/" * label * "_joint-hist.png", fig)

    ζ_max = maximum(ζ)
    fig = Figure(fontsize = 20)
    ax_b = Axis(fig[1, 1][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Buoyancy, }b/\mathrm{m\,s^{-2}}", aspect = 1)
    ax_ζ = Axis(fig[1, 2][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Vertical vorticity, }\zeta/f", aspect = 1)
    hm_b = heatmap!(ax_b, x/1kilometer, y/1kilometer, b);
    hm_ζ = heatmap!(ax_ζ, x/1kilometer, y/1kilometer, ζ/f; colormap = :coolwarm, colorrange = (-ζ_max/1.5f, ζ_max/1.5f));
    Colorbar(fig[1, 1][1, 2], hm_b)
    Colorbar(fig[1, 2][1, 2], hm_ζ)
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    colsize!(fig.layout, 2, Aspect(1, 1.0))
    resize_to_layout!(fig)
    display(fig)
    save("pretty_things/" * label * "_top-pic.png", fig)

end
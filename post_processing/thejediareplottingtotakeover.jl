using Oceananigans, JLD2, Makie, Printf
using Oceananigans.Units
using CairoMakie
using FFTW
using Statistics

function surface_func(file_label, f, input_labels)
    
    filename = "raw_data/" * file_label * "_BI_xy.jld2"
    file = jldopen(filename)
    iterations = parse.(Int, keys(file["timeseries/t"]))
    return [file["timeseries/t/$iter"] for iter in iterations[11:end]], [f([file["timeseries/" * var_lab * "/$iter"] for var_lab in input_labels]) for iter in iterations[11:end]]

end

function surface_function_stats(file_label, f, input_labels)
    # f is a (broadcasting) function on a tuple X, the value of which are determined by input_labels
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

function surface_stats(file_label, var_label)
    f = x -> x[1]
    input_labels = [var_label]
    return surface_function_stats(file_label, f, input_labels)
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
    t, W², ~, ~ = surface_function_stats(label, 𝐮 -> (𝐮[3].^2), ["u", "v", "w"])
    t, ℬ_surf, ~, ~ = surface_function_stats(label, x -> - x[1] .* x[2], ["w", "b"])
    t, ℬ_vol = buoyancy_flux(label)
    ℬ = ℬ_vol
    t, μ_ζ, σ²_ζ, skew_ζ = surface_stats(label, "ζ")
    t, μ_δ, σ²_δ, skew_δ = surface_stats(label, "δ")
    μ_Ri = t .* 0
    ft = t .* 1e-4

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"ft", ylabel=L"\text{Skewness}", width=250, height=250)
    l1 = lines!(ft, skew_ζ, label = L"\zeta")
    l2 = lines!(ft, skew_δ, label = L"\delta")
    axislegend(position = :rb)
    resize_to_layout!(fig)
    save("pretty_things/" * label * "_skew.pdf", fig)
    display(fig)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"ft", ylabel=L"\ln(\sigma^2)", width=250, height=250)
    l1 = lines!(ft, log.(σ²_ζ/f^2), label = L"\zeta/f")
    l2 = lines!(ft, log.(σ²_δ/f^2), label = L"\delta/f")
    l4 = lines!(ft, -log.(σ²_ζ./σ²_δ), label = L"\delta/\zeta")
    axislegend(position = :rt)
    resize_to_layout!(fig)
    save("pretty_things/" * label * "_var.pdf", fig)
    display(fig)
    fig = Figure()
    power_ten = floor(log10(maximum(ℬ)))+6
    @info power_ten
    ax = Axis(fig[1, 1], xlabel = L"ft", ylabel=L"\text{Buoyancy flux, }\mathcal{B}=-\langle w'b'\rangle/\mathrm{mm^2\,s^{-3}}", width=250, height=250)#, yticks = [0, 10^power_ten])
    l = lines!(ft, ℬ/10^(-6))
    resize_to_layout!(fig)
    save("pretty_things/" * label * "_ℬ.pdf", fig)
    display(fig)

end

#=function f_timeseries(label, Ri)    # (Fourier series)

    function abs_f_cpts(field, m_max = 4, l_max = 2)
        fmags = zeros(1+m_max, 1+l_max)
        ffield = fft(field)/length(field)
        M, L = size(ffield)
        for m = 0 : m_max, l = 0 : l_max
            ms = (m == 0 ? [1] : [m+1, M+1-m])
            ls = (l == 0 ? [1] : [l+1, L+1-l])
            fmags[m+1, l+1] = sum((abs.(ffield[ms, ls])).^2)
        end
        return fmags
    end

    f = 1e-4
    s = 1e4
    H = 50
    t, W², ~, ~ = surface_function_stats(label, 𝐮 -> (𝐮[3].^2), ["u", "v", "w"])
    t, fmags = surface_func(label, x -> - abs_f_cpts(x[1]), ["w"])
    for i in length(t)
        for fmag in fmags[i]
            if fmag < 0
                @info i, fmag
            end
        end
    end
    @info typeof(fmags)
    ft = t .* 1e-4

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"ft", ylabel=L"\text{log-RMS}")
    for m = 0:3, l = 0:2
        if m >= 2
            lines!(ft, log.(abs.([fmags[i][m+1, l+1] for i = 1:length(fmags)])) .+ 20, label = string(m)*", "*string(l))
        end
    end
    lines!(ft, log.(W²) .+ 20, label = "total", linestyle = :dash, color = :black)
    axislegend(position = :rc)
    slope = (10.8*(1+Ri))^(-0.5)
    i₁ = Int(round(length(ft)*0.05))
    i₂ = Int(round(length(ft)*0.6))
    lines!(ft[i₁:i₂], ft[i₁:i₂] * 0.71*2*slope .- 30, linestyle = :dot, color = :black)
    save("pretty_things/" * label * "_var.pdf", fig)
    display(fig)

end=#

function get_iter(label, ft)
    f = 1e-4
    filename = "raw_data/" * label * "_BI_xy.jld2"
    file = jldopen(filename)
    iterations = parse.(Int, keys(file["timeseries/t"]))
    ts = [file["timeseries/t/$iter"] for iter in iterations]
    i = argmin(abs.(ts .- ft/f))
    return iterations[i]
end

function ζ_δ_joint_freq(ζ, δ)

    f = 1e-4
    ζ = vec(ζ)
    δ = vec(δ)
    n = 200
    a = -10f
    b = 10f
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

function ζ_δ_joint_hist(label, ft)
    
    f = 1e-4
    iter = get_iter(label, ft)
    filename = "raw_data/" * label * "_BI_xy.jld2"
    file = jldopen(filename)
    ζ = (file["timeseries/ζ/$iter"])
    δ = (file["timeseries/δ/$iter"])
    x, y, freqs = ζ_δ_joint_freq(ζ, δ)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\zeta/f", ylabel = L"\delta/f", width = 250, height = 250, title=L"ft=%$(ft)")
    heatmap!(ax, x/f, y/f, freqs.^0.4, colormap = :amp, colorrange = (0, 20))
    lines!(ax, [0, 0], [-10, 10], color = :gray, linestyle = :dash, linewidth = 1)
    lines!(ax, [-10, 10], [0, 0], color = :gray, linestyle = :dash, linewidth = 1)
    lines!(ax, [-1, -1], [-10, 10], color = :blue, linestyle = :dot, linewidth = 1)
    resize_to_layout!(fig)
    display(fig)

end

function ζ_δ_joint_hists(label, ft1, ft2)
    
    f = 1e-4
    iter1 = get_iter(label, ft1)
    iter2 = get_iter(label, ft2)
    filename = "raw_data/" * label * "_BI_xy.jld2"
    file = jldopen(filename)
    ζ1 = (file["timeseries/ζ/$iter1"])
    δ1 = (file["timeseries/δ/$iter1"])
    ζ2 = (file["timeseries/ζ/$iter2"])
    δ2 = (file["timeseries/δ/$iter2"])
    x, y, freqs1 = ζ_δ_joint_freq(ζ1, δ1)
    ~, ~, freqs2 = ζ_δ_joint_freq(ζ2, δ2)
    fig = Figure()
    ax1 = Axis(fig[1, 1], xlabel = L"\zeta/f", ylabel = L"\delta/f", width = 250, height = 250, title=L"ft=%$(ft1)")
    heatmap!(ax1, x/f, y/f, freqs1.^0.4, colormap = :amp, colorrange = (0, 20))
    lines!(ax1, [0, 0], [-10, 10], color = :gray, linestyle = :dash, linewidth = 1)
    lines!(ax1, [-10, 10], [0, 0], color = :gray, linestyle = :dash, linewidth = 1)
    lines!(ax1, [-1, -1], [-10, 10], color = :brown, linestyle = :dot, linewidth = 1)
    ax2 = Axis(fig[1, 2], xlabel = L"\zeta/f", width = 250, height = 250, title=L"ft=%$(ft2)")
    heatmap!(ax2, x/f, y/f, freqs2.^0.4, colormap = :amp, colorrange = (0, 20))
    lines!(ax2, [0, 0], [-10, 10], color = :gray, linestyle = :dash, linewidth = 1)
    lines!(ax2, [-10, 10], [0, 0], color = :gray, linestyle = :dash, linewidth = 1)
    lines!(ax2, [-1, -1], [-10, 10], color = :brown, linestyle = :dot, linewidth = 1)
    resize_to_layout!(fig)
    display(fig)
    save("pretty_things/" * label * "_joint-hists.png", fig)

end

function N²(label)  # Plot average N² vs. time
    
    f = 1e-4
    filename_y_mean = "raw_data/" * label * "_BI_y-avg.jld2"
    file = jldopen(filename_y_mean)
    iterations = parse.(Int, keys(file["timeseries/t"]))
    t = [file["timeseries/t/$iter"] for iter in iterations]
    avg_N² = [mean(file["timeseries/b̅/$iter"][:, 1, 1] - file["timeseries/b̅/$iter"][:, 1, 64])/50 for iter in iterations]
    lines(f*t, avg_N²)

end
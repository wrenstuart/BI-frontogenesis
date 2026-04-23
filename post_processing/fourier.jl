using FFTW
using Oceananigans
using JLD2
using CairoMakie
include("pp-io.jl")

# if x = 0:0.1:2π and y = x then fft(sin.(x)*sin.(2y)') looks like:
#
#                    0 0 0 0 … 0 0 0
#                    0 0 * 0 … 0 * 0
#                    ...............
#                    0 0 0 0 … 0 0 0
#                    0 0 * 0 … 0 * 0

function k_space(x)
    N = length(x)
    L = (x[end] - x[1]) * N/(N-1)
    k = zeros((N))
    for i = 1 : Int(floor(N/2))
        k[N+1-i] = -i*2π/L
        k[i+1] = i*2π/L
    end
    return k
end

snapshot(file::JLD2.JLDFile{JLD2.MmapIO}, iter, var) = file["timeseries/"*var*"/"*string(iter)][4:end-3,4:end-3,1]
snapshot(label::String, iter, var) = snapshot(jldopen(data_dir(label)*"BI_xy_top.jld2"), iter, var)

function scalar_fft_snapshot(label, ft, var)
    f = 1e-4
    filename = data_dir(label)*"BI_xy_top.jld2"
    x, y, ~ = nodes(FieldTimeSeries(filename, var, iterations = 0))
    file = jldopen(filename)
    iterations = parse.(Int, keys(file["timeseries/t"]))
    ts = [file["timeseries/t/$iter"] for iter in iterations]
    iter = iterations[argmin(abs.(f*ts.-ft))]
    fft_snapshot = fft(snapshot(file, iter, var))
    𝑘_space = k_space(x)
    𝑙_space = k_space(y)
    M, N = size(fft_snapshot)
    k_max = 2^0.5 * maximum(𝑘_space)
    Δk = 2k_max/N
    ks = Array(0:Δk:k_max)
    fs = ks*0
    for i = 1 : M, j = 1 : N
        k = (𝑘_space[i]^2 + 𝑙_space[j]^2) ^ 0.5
        k_index = argmin(abs.(ks .- k))
        fs[k_index] += abs.(fft_snapshot[i, j])^2
    end
    return ks, fs
end

function fft_crosssection_snapshot(label, ft, var)
    f = 1e-4
    filename = data_dir(label)*"BI_xy_top.jld2"
    x, y, ~ = nodes(FieldTimeSeries(filename, var, iterations = 0))
    file = jldopen(filename)
    iterations = parse.(Int, keys(file["timeseries/t"]))
    ts = [file["timeseries/t/$iter"] for iter in iterations]
    iter = iterations[argmin(abs.(f*ts.-ft))]
    crosssection = snapshot(file, iter, var)[:, 1]
    fft_crossection = fft(crosssection)
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, k_space(x), abs.(fft_crossection).^2)
    display(fig)
end
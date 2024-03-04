using Oceananigans, JLD2, Plots, Printf
using Oceananigans.Units

function ζ_var_timeseries(label)
    filename_xy = "raw_data/" * label * "_BI_xy"
    file_xy = jldopen(filename_xy * ".jld2")
    ζ_ic = FieldTimeSeries(filename_xy * ".jld2", "ζ₃", iterations = 0)
    xζ, yζ, zζ = nodes(ζ_ic)
    iterations = parse.(Int, keys(file_xy["timeseries/t"]))
    t_save = zeros(length(iterations))
    ζ_var = zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t = file_xy["timeseries/t/$iter"]
        t_save[i] = t
        if i > 10       # Ignore the first few iterations due to noise
            ζ_xy = file_xy["timeseries/ζ₃/$iter"]
            ζ_var[i] = sum(ζ_xy .* ζ_xy) / Float64(length(xζ) * length(yζ))
        end
    end

    return t_save, ζ_var

end

function ζ_skew_timeseries(label)
    filename_xy = "raw_data/" * label * "_BI_xy"
    file_xy = jldopen(filename_xy * ".jld2")
    ζ_ic = FieldTimeSeries(filename_xy * ".jld2", "ζ₃", iterations = 0)
    xζ, yζ, zζ = nodes(ζ_ic)
    iterations = parse.(Int, keys(file_xy["timeseries/t"]))
    t_save = zeros(length(iterations))
    ζ_skew= zeros(Float64, length(iterations))
    for i = 11 : length(iterations)
        iter = iterations[i]
        t = file_xy["timeseries/t/$iter"]
        t_save[i] = t
        ζ_xy = file_xy["timeseries/ζ₃/$iter"]
        ζ_var = sum(ζ_xy .^ 2) / (length(xζ) * length(yζ))
        ζ_scale = ζ_var ^ 0.5
        ζ_skew[i] = sum(ζ_xy .^ 3) / (length(xζ) * length(yζ) * ζ_scale ^ 3)
    end

    return t_save, ζ_skew

end

f = 1e-4
##################
label = "test1000"
##################
t, ζ_var = ζ_var_timeseries(label)
~, ζ_skew = ζ_skew_timeseries(label)
t = t/t[end] * 40

ζ_max_skew = maximum(ζ_skew)
i_break = findall(x -> x == ζ_max_skew, ζ_skew)[1]

filename_xy = "raw_data/" * label * "_BI_xy"
file_xy = jldopen(filename_xy * ".jld2")
iterations = parse.(Int, keys(file_xy["timeseries/t"]))

iter = iterations[i_break]
ζ_xy = file_xy["timeseries/ζ₃/$iter"][:, :, 1];
ζs = zeros(Float64, (length(ζ_xy) * 5))
ζ_scale = ζ_var[i_break] ^ 0.5
for j = -2 : 2
    @info j
    iter = iterations[minimum([i_break+j, length(iterations)])]
    ζ_xy = file_xy["timeseries/ζ₃/$iter"][:, :, 1];
    ζs[length(ζ_xy) * (j+2) + 1 : length(ζ_xy) * (j+3)] = vec(ζ_xy)
end
@info ζ_skew[i_break]
plt = histogram(ζs/f, nbins = 200, xlims = (-5*ζ_scale/f, 5*ζ_scale/f), label = "", size = (640, 360))
Plots.xlabel!("\$\\zeta/f\$")
title!(label)
display(plt)
savefig("pretty_things/" * label * ".png")
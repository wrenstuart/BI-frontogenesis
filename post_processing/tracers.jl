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
    grid_x, grid_y, grid_z = nodes(FieldTimeSeries(filename_xy_top, "ζ", iterations = 0))
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
    f = (1-x_frac) * (1-y_frac) * f₋₋ + (1-x_frac) * y_frac * f₋₊ + x_frac * (1-y_frac) * f₊₋ + x_frac * y_frac * f₊₊
    
    return f

end

function grid_nearest(data::FileData, var::String, x::Float64, y::Float64, iter::Int)

    i₋ = Int(floor(x/data.Lx * data.Nx)) + 1
    i₊ = i₋ % data.Nx + 1
    j₋ = Int(floor(y/data.Ly * data.Ny)) + 1
    j₊ = j₋ % data.Ny + 1
    if x - data.x[i₋] < data.x[i₊] - x
        i = i₋
    else
        i = i₊
    end
    if y - data.y[j₋] < data.y[j₊] - y
        j = j₋
    else
        j = j₊
    end

    return data.file["timeseries/$var/$iter"][i, j, 1]

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

function lagr_track(data::FileData, var_to_track::Tuple{Function, Vector{String}}, drifter::Vector{Vector{Float64}})
    
    # Tracks the value of the function var_to_track[1] on variables with labels given by var_to_track[2]
    # along the trajectory of drifter
    
    f, input_labels = var_to_track
    iterations = parse.(Int, keys(data.file["timeseries/t"]))
    t = [data.file["timeseries/t/$iter"] for iter in iterations]
    output = zeros(Float64, length(iterations))
    for (i, iter) in enumerate(iterations)
        t[i] = data.file["timeseries/t/$iter"]
        x, y = drifter[i]
        input = [grid_nearest(data, var, x, y, iter) for var in input_labels]
        output[i] = f(input)
    end

    return t, output

end

function lagr_track(data::FileData, var_to_track::String, drifter::Vector{Vector{Float64}})

    # Track the value of the variable with label var_to_track

    return lagr_track(data, (x -> x[1], [var_to_track]), drifter)

end

# Define some plottable quantities

function Ri_func(input)
    b_z, u_z, v_z = input
    return b_z / (u_z^2 + v_z^2)
end

function KE_func(input)
    u, v, w = input
    return (u^2 + v^2 + w^2) / 2
end

function ∇ₕb_func(input)
    b_x, b_y = input
    return (b_x^2 + b_y^2) ^ 0.5
end

function ζ_on_f_func(input)
    f = 1e-4
    ζ = input[1]
    return ζ / f
end
function δ_on_f_func(input)
    f = 1e-4
    δ = input[1]
    return δ / f
end

function F_hor_δ_func(input)
    δ, ζ, u_x, v_x = input
    v_y = δ - u_x
    u_y = v_x - ζ
    return - (u_x^2 + 2v_x*u_y + v_y^2)
end
function F_cor_and_pres_δ_func(input)
    ζ, ζ_g = input
    ζ_ag = ζ - ζ_g
    return f * ζ_ag
end
function F_vert_δ_func(input)
    u_z, v_z, w_x, w_y = input
    return - (u_z * w_x + v_z * w_y)
end
F_hor_δ_appr_func(input) = -input[1]^2
F_cor_and_pres_δ_appr_func(input) = F_cor_and_pres_δ_func(input)
F_vert_δ_appr_func(input) = 0

function F_hor_ζ_func(input)
    δ, ζ = input
    return - δ * ζ
end
function F_cor_ζ_func(input)
    δ = input[1]
    return - f * δ
end
function F_vert_ζ_func(input)
    u_z, v_z, w_x, w_y = input
    return u_z * w_y - v_z * w_x
end
F_hor_ζ_appr_func(input) = F_hor_ζ_func(input)
F_cor_ζ_appr_func(input) = F_cor_ζ_func(input)
F_vert_ζ_appr_func(input) = 0

plotting_vars = (Ri = (Ri_func, ["b_z", "u_z", "v_z"]),
                 KE = (KE_func, ["u", "v", "w"]),
                 ∇ₕb = (∇ₕb_func, ["b_x", "b_y"]),
                 ζ_on_f = (ζ_on_f_func, ["ζ"]),
                 δ_on_f = (δ_on_f_func, ["δ"]),

                 F_hor_δ = (F_hor_δ_func, ["δ", "ζ", "u_x", "v_x"]),
                 F_hor_δ_appr = (F_hor_δ_appr_func, ["δ"]),
                 F_cor_and_pres_δ = (F_cor_and_pres_δ_func, ["ζ", "ζ_g"]),
                 F_cor_and_pres_δ_appr = (F_cor_and_pres_δ_appr_func, ["ζ", "ζ_g"]),
                 F_vert_δ = (F_vert_δ_func, ["u_z", "v_z", "w_x", "w_y"]),
                 F_vert_δ_appr = (F_vert_δ_appr_func, ["δ"]),

                 F_hor_ζ = (F_hor_ζ_func, ["δ", "ζ"]),
                 F_hor_ζ_appr = (F_hor_ζ_appr_func, ["δ", "ζ"]),
                 F_cor_ζ = (F_cor_ζ_func, ["δ"]),
                 F_cor_ζ_appr = (F_cor_ζ_appr_func, ["δ"]),
                 F_vert_ζ = (F_vert_ζ_func, ["u_z", "v_z", "w_x", "w_y"]),
                 F_vert_appr_ζ = (F_vert_ζ_appr_func, ["ζ"]))

function tracer_track(label::String, var_to_track::Union{Tuple{Function, Vector{String}}, String})

    # var_to_track can either be a tuple [f, var_labels], wherein f is evaluated on variables with
    # labels given by var_labels, or it can be a string, in which case just the variable with that
    # label is evaluated

    f = 1e-4

    data = topdata(label)
    drifters = extract_tracers(label)[1:20]
    fig = Figure()
    ax = Axis(fig[1, 1])#, limits = (nothing, (-1, 1)))
    for drifter in drifters
        t, var = lagr_track(data, var_to_track, drifter)
        interesting_bits(data, drifter)
        #=i₁ = Int(round(length(t)/3))
        i₂ = Int(round(2length(t)/3))=#
        i₁ = 5
        i₂ = length(t) - 4
        var_smooth = zeros(length(var))
        for i = i₁ : i₂
            var_smooth[i] = sum(var[i-4:i+4])/9
        end
        lines!(ax, f*t[i₁:i₂], var_smooth[i₁:i₂])
        #lines!(ax, t[i₁:i₂], var[i₁:i₂])
    end
    display(fig)
    #save("pretty_things/tracer-delta_" * label * ".pdf", fig)
    
end

function interesting_bits(data::FileData, drifter)
    
    f = 1e-4
    t, ζ = lagr_track(data, "ζ", drifter)
    t, δ = lagr_track(data, "δ", drifter)
    t, u_x = lagr_track(data, "u_x", drifter)
    t, v_x = lagr_track(data, "v_x", drifter)
    i₁ = 5
    i₂ = length(t) - 4
    t = t[i₁:i₂]
    ζ = sum([ζ[i₁+j : i₂+j] for j = -4 : 4])
    δ = sum([δ[i₁+j : i₂+j] for j = -4 : 4])
    u_x = sum([u_x[i₁+j : i₂+j] for j = -4 : 4])
    v_x = sum([v_x[i₁+j : i₂+j] for j = -4 : 4])
    u_y = v_x - ζ
    v_y = δ - u_x
    abs_∇𝐮 = (u_x.^2 + u_y.^2 + v_x.^2 + v_y.^2) .^ 0.5
    sections = []
    run_start = 1
    for i = 1 : length(t)
        if abs_∇𝐮[i] > 5f
            if run_start == i
                push!(sections, [])
            end
            push!(sections[end], i)
        else
            run_start = i+1
        end
    end

    #=if length(sections) > 0

        s = sections[1]
        t₀ = t[argmax(abs.(ζ[s])) + s[1] - 1]
        if f*t₀ < 20 && length(sections) > 1
            s = sections[2]
            t₀ = t[argmax(abs.(ζ[s])) + s[1] - 1]
        end
        T = minimum([t[s[end]] - t₀, t₀ - t[s[1]]])
        t_rel = (t[s].-t₀)/T
        
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel = L"\zeta/f", ylabel = L"\delta/f")
        lines!(ax, ζ[s]/f, δ[s]/f, color = t_rel, colormap = :coolwarm, colorrange = (-1.0,1.0))
        
        #ax = Axis(fig[1, 1], xlabel = L"ft", ylabel = L"\zeta/f,\:\delta/f")
        #lines!(ax, f*(t[s1].-t₀), ζ[s1]/f)
        #lines!(ax, f*(t[s1].-t₀), δ[s1]/f)

        display(fig)
    end=#

    trajectories = []
    for s in sections
        t₀ = t[argmax(abs.(ζ[s])) + s[1] - 1]
        if f*t₀ < 20 && length(sections) > 1
            s = sections[2]
            t₀ = t[argmax(abs.(ζ[s])) + s[1] - 1]
        end
        T = minimum([t[s[end]] - t₀, t₀ - t[s[1]]])
        t_rel = (t[s].-t₀)/T
        push!(trajectories, (ζ = ζ[s], δ = δ[s], t_rel = t_rel))
    end

    return trajectories

end

function all_traj_ζ_δ_plot(label::String)
    
    f = 1e-4

    data = topdata(label)
    drifters = extract_tracers(label)
    trajectories = []
    for drifter in drifters
        [push!(trajectories, traj) for traj in interesting_bits(data, drifter)]
    end

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\zeta/f", ylabel = L"\delta/f")
    [lines!(ax, traj.ζ/f, traj.δ/f, color = traj.t_rel, colormap = :coolwarm, colorrange = (-1.0,1.0)) for traj in trajectories]
    display(fig)

end
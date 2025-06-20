using CairoMakie
using Oceananigans
using Oceananigans.Units
using JLD2
using Printf
using Colors
using Statistics

include("tracers.jl")

δ_mult = 1
other_mult = 1
νₕ = 1e+1
νᵥ = 1e-3

function rearrange_data(label)
    in_name = "raw_data/" * label * "_particles.jld2"
    out_name = "raw_data/" * label * "_particles-rearranged.jld2"
    if ~isfile(in_name)
        throw("Could not find drifter data associated with label \'" * label * "\'")
    end
    in_file = jldopen(in_name)
    iterations = parse.(Int, keys(in_file["timeseries/t"]))
    t = [in_file["timeseries/t/$iter"] for iter in iterations]
    raw_keys = keys(in_file["timeseries/particles/0"])
    N = length(in_file["timeseries/particles/0"].x)
        # Number of particles
    p_data = in_file["timeseries/particles"]
    drifters = [[(; (raw_keys .=> [p_data["$iter"][var][n] for var in raw_keys])...) for iter in iterations] for n = 1 : N]
    @save out_name t drifters
end

function drifter_keys(label)
    filename = "raw_data/" * label * "_particles.jld2"
    file = jldopen(filename)
    return keys(file["timeseries/particles/0"])
end

function get_drifter_data(label)
    filename = "raw_data/" * label * "_particles-rearranged.jld2"
    if ~isfile(filename)
        @info "Rearranging output data..."
        rearrange_data(label)
    end
    file = jldopen(filename)
    t = file["t"]
    drifters = file["drifters"]
    return t, drifters
end

label = "test"
t, drifters = get_drifter_data(label)

f = 1e-4

# In the following functions, d is the data associated
# with a single particle at a given iteration
u_y(d) = d.v_x - d.ζ
v_y(d) = d.δ - d.u_x
∇ₕb(d) = (d.b_x^2 + d.b_y^2) ^ 0.5
F_hor_ζ(d) = d.F_ζ_hor#-d.δ * d.ζ
F_vrt_ζ(d) = d.F_ζ_vrt#d.w_y * d.u_z - d.w_x * d.v_z
F_Cor_ζ(d) = d.ζ_cor#-f * d.δ
H_mix_ζ(d) = d.ζ_visc#νₕ * d.∇ₕ²ζ
V_mix_ζ(d) = 0#νᵥ * d.ζ_zz
ζ_err_func(d) = d.ζ_err
vert_adv_ζ(d) = 0#-d.w * d.ζ_z
ζ_adv_func(d) = d.ζ_adv
ζ_tendency_func(d) = d.ζ_tendency

F_hor_δ(d) = 0#-(d.u_x ^ 2 + 2 * d.v_x * u_y(d) + v_y(d) ^ 2)
F_vrt_δ(d) = 0#-(d.w_x * d.u_z + d.w_y * d.v_z)
F_Cor_δ(d) = 0#f * d.ζ
F_prs_δ(d) = 0#-d.fζ_g
H_mix_δ(d) = 0#νₕ * d.∇ₕ²δ
V_mix_δ(d) = 0#νᵥ * d.δ_zz

∇ₕ𝐮ₕ(d) = (d.u_x ^ 2 + d.v_x ^ 2 + u_y(d) ^ 2 + v_y(d) ^ 2) ^ 0.5
ζ_on_f(d) = d.ζ / f
δ_on_f(d) = d.δ / f

function get_sections(drifter)
    abs_∇ₕ𝐮ₕ = ∇ₕ𝐮ₕ.(drifter)
    sections = []
    run_start = 2
    for i = 2 : length(t)-1
        if abs_∇ₕ𝐮ₕ[i] > 5f
            if run_start == i
                push!(sections, [])
            end
            push!(sections[end], i)
        else
            run_start = i+1
        end
    end
    return sections
end

function ζ_δ_trajectories(drifter)
    # Finds segements of a drifter's trajectory through ζ–δ phase space that are high in ∇𝐮 and
    # return an array of each trajectory, where each trajectory is a named tuple with keys ζ, δ and t_rel
    # t_rel is a normalised time relative to the time of maximal ζ

    ζ_f = ζ_on_f.(drifter)
    δ_f = δ_on_f.(drifter)
    sections = get_sections(drifter)

    #=if length(sections) > 0

        s = sections[1]
        t₀ = t[argmax(abs.(ζ_f[s])) + s[1] - 1]
        if f*t₀ < 20 && length(sections) > 1
            s = sections[2]
            t₀ = t[argmax(abs.(ζ_f[s])) + s[1] - 1]
        end
        T = minimum([t[s[end]] - t₀, t₀ - t[s[1]]])
        t_rel = (t[s].-t₀)/T
        
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel = L"\zeta/f", ylabel = L"\delta/f")
        lines!(ax, ζ_f[s], δ_f[s], color = t_rel, colormap = :coolwarm, colorrange = (-1.0,1.0))

        #fig = Figure()
        #ax = Axis(fig[1, 1], xlabel = L"ft")
        #lines!(ax, f*t[s], δ[s])
        #lines!(ax, f*t[s], ζ[s])

        display(fig)

    end=#

    trajectories = []
    for s in sections
        t₀ = t[argmax(abs.(ζ_f[s])) + s[1] - 1]
        if f*t₀ < 20 && length(sections) > 1
            s = sections[2]
            t₀ = t[argmax(abs.(ζ_f[s])) + s[1] - 1]
        end
        T = minimum([t[s[end]] - t₀, t₀ - t[s[1]]])
        t_rel = (t[s].-t₀)/T
        push!(trajectories, (ζ = ζ_f[s], δ = δ_f[s], t_rel = t_rel))
    end

    return trajectories

end

function ζ_δ_arrow_map(drifters, ζ_grid = -2:0.5:20, δ_grid = -10:0.5:5)#ζ_grid = -5:0.5:49, δ_grid = -49:0.5:5)

    n_ζ = length(ζ_grid)
    n_δ = length(δ_grid)
    ζ₁ = ζ_grid[1]
    ζₙ = ζ_grid[end]
    δ₁ = δ_grid[1]
    δₙ = δ_grid[end]
    ζ_lims = (2ζ₁-ζ_grid[2], 2ζₙ-ζ_grid[end-1])
    δ_lims = (2δ₁-δ_grid[2], 2δₙ-δ_grid[end-1])

    total_sections_count = 0
    visited_count_total = zeros((n_ζ, n_δ))
    time_count_total = zeros((n_ζ, n_δ))
    Dₜζ_total = zeros((n_ζ, n_δ))
    Dₜδ_total = zeros((n_ζ, n_δ))

    for drifter in drifters

        #for sec in get_sections(drifter)
        for sec in [2 : length(t)-1]

            total_sections_count += 1
            visited_count = zeros(Int64, (n_ζ, n_δ))
            ζ_f = ζ_on_f.(drifter[sec])
            δ_f = δ_on_f.(drifter[sec])
            Dₜζ_f = (ζ_on_f.(drifter[sec.+1]) - ζ_on_f.(drifter[sec.-1])) ./ (t[sec.+1] - t[sec.-1])
            Dₜδ_f = (δ_on_f.(drifter[sec.+1]) - δ_on_f.(drifter[sec.-1])) ./ (t[sec.+1] - t[sec.-1])
            for (j, ζ) in enumerate(ζ_f)
                δ = δ_f[j]
                if ζ_lims[1] < ζ && ζ < ζ_lims[2] && δ_lims[1] < δ && δ < δ_lims[2]
                    i₁ = argmin(abs.(ζ .- ζ_grid))
                    i₂ = argmin(abs.(δ .- δ_grid))
                    visited_count[i₁, i₂] |= 1
                    time_count_total[i₁, i₂] += 1
                    Dₜζ_total[i₁, i₂] += Dₜζ_f[j]
                    Dₜδ_total[i₁, i₂] += Dₜδ_f[j]
                end
            end
            visited_count_total += visited_count

        end

    end

    avg_weight = [n == 0 ? 0 : 1/n for n in time_count_total]
    Dₜζ̅ = Dₜζ_total .* avg_weight# * [n > 10 ? 1 : 0 for n in visited_count_total]
    Dₜδ̅ = Dₜδ_total .* avg_weight# * [n > 10 ? 1 : 0 for n in visited_count_total]
    #arrows(ζ_grid, δ_grid, Dₜζ̅, Dₜδ̅, lengthscale = 0.1)
    mags = (Dₜζ̅.^2 + Dₜδ̅.^2) .^ 0.5
    mags = [x == 0 ? 1 : x for x in mags]
    #arrows(ζ_grid, δ_grid, Dₜζ̅, Dₜδ̅, lengthscale = 100)

    mags = (Dₜζ_total.^2 + Dₜδ_total.^2) .^ 0.5
    u = Dₜζ_total./mags
    v = Dₜδ_total./mags
    u .*= mags.^0.2
    v .*= mags.^0.2
    #alpha = log.(1 .+ visited_count_total)
    alpha = visited_count_total#.^0.5
    colors = [RGBA(0, 0, 0, α/maximum(alpha)) for α in alpha]
    arrows(ζ_grid, δ_grid, u, v,
           lengthscale = 0.35, color = vec(colors), lengths = mags)

end

function lagr_ζ_balance(drifter)
    sections = get_sections(drifter)
    for s in sections
        if length(s) > 1
            plot_ζ_balance(drifter, s)
        end
    end
end

function lagr_δ_balance(drifter)
    sections = get_sections(drifter)
    for s in sections
        if length(s) > 1
            plot_δ_balance(drifter, s)
        end
    end
end

function ζ_δ(drifter)
    sections = get_sections(drifter)
    for s in sections
        if length(s) > 1
            ζ_δ(drifter, s)
        end
    end
end

function smooth_timeseries(x)
    if length(x) < 5
        return x
    else
        return [mean(x[i-2:i+2]) for i = 3:length(x)-2]
    end
end

function plot_ζ_balance(drifter, section)
    s = section
    Δt = t[s[2]] - t[s[1]]
    ft = f * t[s]
    Dₜζ = [(drifter[i+1].ζ - drifter[i-1].ζ) / (2Δt) for i in s]*other_mult
    fig = Figure()
    ax = Axis(fig[1, 1])
    #ylims!(ax, (-2e-7, 2e-7))
    F_hor = F_hor_ζ.(drifter[s])*δ_mult*other_mult
    F_vrt = F_vrt_ζ.(drifter[s])*δ_mult*other_mult
    F_Cor = F_Cor_ζ.(drifter[s])*δ_mult
    H_mix = H_mix_ζ.(drifter[s])*other_mult
    V_mix = V_mix_ζ.(drifter[s])*other_mult
    vert_adv = vert_adv_ζ.(drifter[s])
    if length(ft) > 5
        ft = ft[3:end-2]
        Dₜζ = smooth_timeseries(Dₜζ)
        F_hor = smooth_timeseries(F_hor)
        F_vrt = smooth_timeseries(F_vrt)
        F_Cor = smooth_timeseries(F_Cor)
        H_mix = smooth_timeseries(H_mix)
        V_mix = smooth_timeseries(V_mix)
        vert_adv = smooth_timeseries(vert_adv)
    end
    lines!(ax, ft, F_hor, label = "horizontal")
    lines!(ax, ft, F_Cor, label = "Coriolis")
    lines!(ax, ft, H_mix, label = "hor. mixing")
    lines!(ax, ft, V_mix, label = "vert. mixing")
    lines!(ax, ft, vert_adv, label = "vert. adv.")
    lines!(ax, ft, Dₜζ - (F_hor + F_vrt + F_Cor + H_mix + V_mix + vert_adv), label = "discrepancy", color = :black, linestyle = :dot)
    lines!(ax, ft, Dₜζ, label = L"D\zeta/Dt", color = :black)
    axislegend()
    display(fig)
end

function ζ_δ(drifter, section)
    s = section
    Δt = t[s[2]] - t[s[1]]
    ft = f * t[s]
    Dₜζ = [(drifter[i+1].ζ - drifter[i-1].ζ) / (2Δt) for i in s]*other_mult
    fig = Figure()
    ax = Axis(fig[1, 1])
    #ylims!(ax, (-2e-7, 2e-7))
    ζ = ζ_on_f.(drifter[s])
    δ = δ_on_f.(drifter[s])
    ζ_g = -F_prs_δ.(drifter[s])/f
    if length(ft) > 5
        ft = ft[3:end-2]
        ζ = smooth_timeseries(ζ)
        δ = smooth_timeseries(δ)
        ζ_g = smooth_timeseries(ζ_g)
    end
    lines!(ax, ft, ζ, label = "zeta")
    lines!(ax, ft, δ, label = "delta")
    lines!(ax, ft, 1000ζ_g, label = "zeta_g")
    axislegend()
    display(fig)
end

#=function investigate(drifter)

    ft = f*t
    fig = Figure()
    ax = Axis(fig[1, 1], limits = ((10, 30), (-0.0005, 0.0005)))
    δ_app = [d.δ for d in drifter]
    ζ_zz_app = [d.ζ_zz for d in drifter]
    ζ = [d.ζ for d in drifter]
    ζ_zz_true = 2*(ζ_zz_app+0.5*ζ/Δz^2)
    x = [d.y for d in drifter]

    # now for the interpolation bit

    data = topdata(label)
    drifter = extract_tracers(label)[1]
    new_t, new_δ = lagr_track_new(label, data, "δ", drifter)

    lines!(ax, ft, δ_app)
    lines!(ax, f*new_t, 0.5*new_δ)
    display(fig)

end=#

function plot_δ_balance(drifter, section)       # needs fixing!!!! iterations and t don't line up
    s = section
    Δt = t[s[2]] - t[s[1]]
    ft = f * t[s]
    Dₜδ = [(drifter[i+1].δ - drifter[i-1].δ) / (2Δt) for i in s]*δ_mult
    fig = Figure()
    ax = Axis(fig[1, 1])
    F_hor = F_hor_δ.(drifter[s])*δ_mult^2
    F_Cor = F_Cor_δ.(drifter[s])*other_mult
    F_prs = F_prs_δ.(drifter[s])*other_mult
    H_mix = H_mix_δ.(drifter[s])*δ_mult
    V_mix = V_mix_δ.(drifter[s])*δ_mult
    if length(ft) > 5
        ft = ft[3:end-2]
        Dₜδ = smooth_timeseries(Dₜδ)
        F_hor = smooth_timeseries(F_hor)
        F_Cor = smooth_timeseries(F_Cor)
        F_prs = smooth_timeseries(F_prs)
        H_mix = smooth_timeseries(H_mix)
        V_mix = smooth_timeseries(V_mix)
    end
    lines!(ax, ft, F_hor, label = "horizontal")
    lines!(ax, ft, F_Cor, label = "Coriolis")
    lines!(ax, ft, F_prs, label = "pressure")
    lines!(ax, ft, H_mix, label = "hor. mixing")
    lines!(ax, ft, V_mix, label = "vert. mixing")
    lines!(ax, ft, Dₜδ - (F_hor + F_Cor + H_mix + V_mix + F_prs), label = "discrepancy", color = :black, linestyle = :dot)
    lines!(ax, ft, Dₜδ, label = L"D\delta/Dt", color = :black)
    axislegend(position = :lb)
    display(fig)
end

function animate_drifter_balance(drifter, section)

    fig = Figure(size = (950, 950))
    frame = Observable(1)

    if length(section) < 5
        return
    end
    s = section
    Δt = t[s[2]] - t[s[1]]
    ft = f * t[s]
    Dₜζ = [(drifter[i+1].ζ - drifter[i-1].ζ) / (2Δt) for i in s]*other_mult
    Dₜδ = [(drifter[i+1].δ - drifter[i-1].δ) / (2Δt) for i in s]*δ_mult
    x = [drifter[i].x for i in s]
    y = [drifter[i].y for i in s]
    Fζ_hor = F_hor_ζ.(drifter[s])*δ_mult*other_mult
    Fζ_vrt = 0 * F_hor_ζ.(drifter[s])          # Due to no-penetration condition
    Fζ_Cor = F_Cor_ζ.(drifter[s])*δ_mult
    Hζ_mix = H_mix_ζ.(drifter[s])*other_mult
    Vζ_mix = V_mix_ζ.(drifter[s])*other_mult
    Fδ_hor = F_hor_δ.(drifter[s])*δ_mult^2
    Fδ_Cor = F_Cor_δ.(drifter[s])*other_mult
    Fδ_prs = F_prs_δ.(drifter[s])*other_mult
    Hδ_mix = H_mix_δ.(drifter[s])*δ_mult
    Vδ_mix = V_mix_δ.(drifter[s])*δ_mult

    s = s[3:end-2]
    ft = ft[3:end-2]
    Dₜζ = smooth_timeseries(Dₜζ)
    Dₜδ = smooth_timeseries(Dₜδ)
    x = smooth_timeseries(x)
    y = smooth_timeseries(y)
    Fζ_hor = smooth_timeseries(Fζ_hor)
    Fζ_Cor = smooth_timeseries(Fζ_Cor)
    Hζ_mix = smooth_timeseries(Hζ_mix)
    Vζ_mix = smooth_timeseries(Vζ_mix)
    Fδ_hor = smooth_timeseries(Fδ_hor)
    Fδ_Cor = smooth_timeseries(Fδ_Cor)
    Fδ_prs = smooth_timeseries(Fδ_prs)
    Hδ_mix = smooth_timeseries(Hδ_mix)
    Vδ_mix = smooth_timeseries(Vδ_mix)

    ft_obs = lift(i -> ft[i], frame)
    x_obs = lift(i -> x[i], frame)
    y_obs = lift(i -> y[i], frame)

    fig = Figure()
    ax_ζ = Axis(fig[2, 1:3], height = 200)
    lines!(ax_ζ, ft, Fζ_hor, label = "horizontal")
    lines!(ax_ζ, ft, Fζ_Cor, label = "Coriolis")
    lines!(ax_ζ, ft, Hζ_mix, label = "hor. mixing")
    lines!(ax_ζ, ft, Vζ_mix, label = "vert. mixing")
    lines!(ax_ζ, ft, Dₜζ - (Fζ_hor + Fζ_Cor + Hζ_mix + Vζ_mix), label = "discrepancy", color = :black, linestyle = :dot)
    lines!(ax_ζ, ft, Dₜζ, label = L"D\zeta/Dt", color = :black)
    vlines!(ax_ζ, ft_obs, color = :black)
    axislegend()

    ax_δ = Axis(fig[3, 1:3], height = 200)
    lines!(ax_δ, ft, Fδ_hor, label = "horizontal")
    lines!(ax_δ, ft, Fδ_Cor, label = "Coriolis")
    lines!(ax_δ, ft, Fδ_prs, label = "pressure")
    lines!(ax_δ, ft, Hδ_mix, label = "hor. mixing")
    lines!(ax_δ, ft, Vδ_mix, label = "vert. mixing")
    lines!(ax_δ, ft, Dₜδ - (Fδ_hor + Fδ_Cor + Hδ_mix + Vδ_mix + Fδ_prs), label = "discrepancy", color = :black, linestyle = :dot)
    lines!(ax_δ, ft, Dₜδ, label = L"D\delta/Dt", color = :black)
    vlines!(ax_δ, ft_obs, color = :black)
    axislegend(position = :lb)

    display(fig)




    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"

    # Read in the first iteration. We do this to load the grid
    b_ic = FieldTimeSeries(filename_xy_top, "b", iterations = 0)
    ζ_ic = FieldTimeSeries(filename_xy_top, "ζ", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xb, yb, ~ = nodes(b_ic)
    xζ, yζ, ~ = nodes(ζ_ic)

    # Now, open the file with our data
    file = jldopen(filename_xy_top)
    # Extract the values that iter can take
    iterations = parse.(Int, keys(file["timeseries/t"]))

    # Set up observables for plotting that will update as the iteration number is updated
    iter = lift(i -> iterations[s[i]], frame)   # Timestep iteration
    ζ_xy = lift(iter -> file["timeseries/ζ/$iter"][:, :, 1], iter) # Surface vertical vorticity at this iteration
    ζ_on_f = lift(iter -> ζ_xy[]/f, iter)
    δ = lift(iter -> file["timeseries/δ/$iter"][:, :, 1], iter)
    δ_on_f = lift(iter -> file["timeseries/δ/$iter"][:, :, 1]/f, iter)
    b_xy = lift(iter -> file["timeseries/b/$iter"][:, :, 1], iter)   # Surface buoyancy at this iteration

    # Calculate the maximum relative vorticity and buoyancy flux to set the scale for the colourmap
    ζ_max = 0
    b_max = maximum(b_ic)

    for i = Int(round(length(iterations)/10)) : length(iterations)
        iter[] = iterations[i]
        ζ_max = maximum([ζ_max, maximum(ζ_xy[])])
    end

    ζ_max = minimum([ζ_max, 20f])

    @info "Drawing first frame"

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    iter[] = 0
    ax_b = Axis(fig[1, 1][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Buoyancy, }b", width = 200, height = 200)
    ax_ζ = Axis(fig[1, 2][1, 1], xlabel = L"$x/\mathrm{km}$", title = L"\text{Vertical vorticity, }\zeta/f", width = 200, height = 200)
    ax_δ = Axis(fig[1, 3][1, 1], xlabel = L"$x/\mathrm{km}$", title = L"\text{Horizontal divergence, }\delta/f", width = 200, height = 200)
    hm_b = heatmap!(ax_b, xb/1kilometer, yb/1kilometer, b_xy; colorrange = (-0.5*b_max, 1.5*b_max));
    #scatter!(ax_b, x_obs, y_obs, marker = '.', markersize = 30, color = :black)
    hm_ζ = heatmap!(ax_ζ, xζ/1kilometer, yζ/1kilometer, ζ_on_f; colormap = :coolwarm, colorrange = (-ζ_max/f, ζ_max/f));
    #scatter!(ax_ζ, x_obs, y_obs, marker = '.', markersize = 30, color = :black, update_limits = false)
    hm_δ = heatmap!(ax_δ, xζ/1kilometer, yζ/1kilometer, δ_on_f; colormap = :coolwarm, colorrange = (-ζ_max/f, ζ_max/f));
    #scatter!(ax_δ, x_obs, y_obs, marker = '.', markersize = 30, color = :black)
    Colorbar(fig[1, 1][1, 2], hm_b, height = 200)
    Colorbar(fig[1, 2][1, 2], hm_ζ, height = 200)
    Colorbar(fig[1, 3][1, 2], hm_δ, height = 200)
    resize_to_layout!(fig)



    display(fig)

    @info "Making an animation from saved data..."
    CairoMakie.record(i -> frame[] = i,
        fig,
        "pretty_things/" * label * ".mp4",
        1 : length(s),
        framerate = 20)

end

function ani_drifters(label::String)
    t_drifters, drifters = get_drifter_data(label)
    ani_drifters(label, drifters[1], 2 : length(t_drifters)-1)
end

function ani_drifters(label::String, drifter, section)     # Animate drifters at the surface over a ζ video
    
    fig = Figure(size = (950, 950))

    s = section
    Δt = t[s[2]] - t[s[1]]
    data = topdata(label)
    iterations = parse.(Int, keys(data.file["timeseries/t"]))
    tracers = extract_tracers(label)
    t_drifters, ~ = get_drifter_data(label)
    ft = f * t_drifters
    first_iter_index = argmin(abs.([data.file["timeseries/t/$iter"] for iter in iterations] .- t_drifters[s[1]])) + 3
    last_iter_index = argmin(abs.([data.file["timeseries/t/$iter"] for iter in iterations] .- t_drifters[s[end]])) - 3

    frame = Observable(1)
    iter = lift(i -> iterations[i], frame)
    t_obs = lift(iter -> data.file["timeseries/t/$iter"], iter)
    ζ_on_f = lift(i -> data.file["timeseries/ζ/$i"][:, :, 1]/f, iter)
    δ_on_f = lift(i -> data.file["timeseries/δ/$i"][:, :, 1]/f, iter)
    b = lift(i -> data.file["timeseries/b/$i"][:, :, 1], iter)
    i_drifter = lift(t -> argmin(abs.(t_drifters .- t)), t_obs)
    ft_obs = lift(i -> ft[i], i_drifter)
    tracers_now_x = lift(i -> drifter[i].x/1e3, i_drifter)
    tracers_now_y = lift(i -> drifter[i].y/1e3, i_drifter)

    if length(section) < 5
        return
    end
    s = section
    ft = f * t[s]
    Dₜζ = [(drifter[i+1].ζ - drifter[i-1].ζ) / (t[i+1] - t[i-1]) for i in s]*other_mult
    Dₜδ = [(drifter[i+1].δ - drifter[i-1].δ) / (t[i+1] - t[i-1]) for i in s]*δ_mult
    x = [drifter[i].x for i in s]
    y = [drifter[i].y for i in s]
    Fζ_hor = F_hor_ζ.(drifter[s])*δ_mult*other_mult
    Fζ_vrt = F_vrt_ζ.(drifter[s])*δ_mult*other_mult
    Fζ_Cor = F_Cor_ζ.(drifter[s])*δ_mult
    Hζ_mix = H_mix_ζ.(drifter[s])*other_mult
    Vζ_mix = V_mix_ζ.(drifter[s])*other_mult
    # new ↓
    ζ_err = ζ_err_func.(drifter[s])
    ζ_adv = ζ_adv_func.(drifter[s])
    ζ_tendency = ζ_tendency_func.(drifter[s])
    # new ↑
    ζ_vert_adv = vert_adv_ζ.(drifter[s])
    Fδ_hor = F_hor_δ.(drifter[s])*δ_mult^2
    Fδ_Cor = F_Cor_δ.(drifter[s])*other_mult
    Fδ_prs = F_prs_δ.(drifter[s])*other_mult
    Hδ_mix = H_mix_δ.(drifter[s])*δ_mult
    Vδ_mix = V_mix_δ.(drifter[s])*δ_mult

    s = s[3:end-2]
    ft = ft[3:end-2]
    Dₜζ = smooth_timeseries(Dₜζ)
    Dₜδ = smooth_timeseries(Dₜδ)
    x = smooth_timeseries(x)
    y = smooth_timeseries(y)
    Fζ_hor = smooth_timeseries(Fζ_hor)
    Fζ_vrt = smooth_timeseries(Fζ_vrt)
    Fζ_Cor = smooth_timeseries(Fζ_Cor)
    Hζ_mix = smooth_timeseries(Hζ_mix)
    Vζ_mix = smooth_timeseries(Vζ_mix)
    ζ_vert_adv = smooth_timeseries(ζ_vert_adv)
    # new ↓
    ζ_err = smooth_timeseries(ζ_err)
    ζ_adv = smooth_timeseries(ζ_adv)
    ζ_tendency = smooth_timeseries(ζ_tendency)
    # new ↑
    Fδ_hor = smooth_timeseries(Fδ_hor)
    Fδ_Cor = smooth_timeseries(Fδ_Cor)
    Fδ_prs = smooth_timeseries(Fδ_prs)
    Hδ_mix = smooth_timeseries(Hδ_mix)
    Vδ_mix = smooth_timeseries(Vδ_mix)

    ax_ζ = Axis(fig[2, 1:3], height = 200)
    lines!(ax_ζ, ft, Fζ_hor, label = "horizontal")
    lines!(ax_ζ, ft, Fζ_vrt, label = "vertical")
    lines!(ax_ζ, ft, Fζ_Cor, label = "Coriolis")
    lines!(ax_ζ, ft, Hζ_mix, label = "mixing")
    lines!(ax_ζ, ft, ζ_err, label = "error")
    lines!(ax_ζ, ft, ζ_vert_adv, label = "vert. adv.")
    # lines!(ax_ζ, ft, Dₜζ - (Fζ_hor + Fζ_vrt + Fζ_Cor + Hζ_mix + Vζ_mix + ζ_vert_adv), label = "discrepancy", color = :black, linestyle = :dot)
    # lines!(ax_ζ, ft, Dₜζ - (Fζ_hor + Fζ_vrt + Fζ_Cor + Hζ_mix + Vζ_mix + ζ_vert_adv + ζ_err), label = "discrepancy", color = :black, linestyle = :dot)
    lines!(ax_ζ, ft, ζ_tendency + ζ_adv - (Fζ_hor + Fζ_vrt + Fζ_Cor + Hζ_mix + Vζ_mix + ζ_vert_adv - ζ_err), label = "discrepancy", color = :black, linestyle = :dot)
    # lines!(ax_ζ, ft, Dₜζ, label = L"D\zeta/Dt", color = :black)
    lines!(ax_ζ, ft, ζ_tendency + ζ_adv, label = L"D\zeta/Dt", color = :black)
    vlines!(ax_ζ, ft_obs, color = :black)
    axislegend(position = :lb)
    # DELTA STUFF HAS NOT BEEN UPDATED
    ax_δ = Axis(fig[3, 1:3], height = 200)
    lines!(ax_δ, ft, Fδ_hor, label = "horizontal")
    lines!(ax_δ, ft, Fδ_Cor, label = "Coriolis")
    lines!(ax_δ, ft, Fδ_prs, label = "pressure")
    lines!(ax_δ, ft, Hδ_mix, label = "hor. mixing")
    lines!(ax_δ, ft, Vδ_mix, label = "vert. mixing")
    lines!(ax_δ, ft, Dₜδ - (Fδ_hor + Fδ_Cor + Hδ_mix + Vδ_mix + Fδ_prs), label = "discrepancy", color = :black, linestyle = :dot)
    lines!(ax_δ, ft, Dₜδ, label = L"D\delta/Dt", color = :black)
    vlines!(ax_δ, ft_obs, color = :black)
    axislegend(position = :lb)

    b_ic = data.file["timeseries/b/0"][:, :, 1]
    b_max = maximum(b_ic)

    ax_ζ = Axis(fig[1, 1][1, 1], aspect = 1)
    hm_ζ = heatmap!(ax_ζ, data.x/1e3, data.y/1e3, ζ_on_f, colormap = :coolwarm, colorrange = (-20, 20), height = 200);
    scatter!(ax_ζ, tracers_now_x, tracers_now_y, marker = '.', markersize = 15, color = :black)
    Colorbar(fig[1, 1][1, 2], hm_ζ, height = 200)
    ax_δ = Axis(fig[1, 2][1, 1], aspect = 1)
    hm_δ = heatmap!(ax_δ, data.x/1e3, data.y/1e3, δ_on_f, colormap = :coolwarm, colorrange = (-20, 20), height = 200);
    scatter!(ax_δ, tracers_now_x, tracers_now_y, marker = '.', markersize = 15, color = :black)
    Colorbar(fig[1, 2][1, 2], hm_δ, height = 200)
    ax_b = Axis(fig[1, 3][1, 1], aspect = 1)
    hm_b = heatmap!(ax_b, data.x/1e3, data.y/1e3, b, colorrange = (-0.5b_max, 1.5b_max), height = 200);
    scatter!(ax_b, tracers_now_x, tracers_now_y, marker = '.', markersize = 15, color = :black)
    Colorbar(fig[1, 3][1, 2], hm_b, height = 200)
    
    resize_to_layout!(fig)
    display(fig)

    record(i -> frame[] = i, fig, "pretty_things/tracer_" * label * ".mp4", first_iter_index : last_iter_index, framerate = 20)
    
end

#ζ_δ_arrow_map(drifters)
#ani_drifters("nu1e2", drifters[24], 2:100)
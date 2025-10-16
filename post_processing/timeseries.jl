include("pp-io.jl")
include("drifters-refactored.jl")

label = "new_zeta_balance2"

#=function investigate_lagr_ζ_balance(label::String, drifter_num::Int64, plot_mode = "tracked")
    
    # WE ASSUME THAT BI_XY AND PARTICLE ITERATIONS ARE THE SAME
    # This plots the various terms affecting the Lagrangian change in ζ,
    # following particles as we do so. "tracked" mode uses particle-tracked data,
    # while "interpolated" looks at the Eulerian fields saved at the top and
    # interpolates them onto the particles in post-processing

    check_pp_dir(label)
    eul_data = topdata(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    num_iters = length(tracked_drifter_data[drifter_num])
    iterations = eul_data.iterations

    grid_pos = (Face(), Face())

    x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]
    y = [tracked_drifter_data[drifter_num][i].y for i = 1 : num_iters]
    ζ_tracked = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]
    ζ_tendency_tracked = [tracked_drifter_data[drifter_num][i].ζ_tendency for i = 1 : num_iters]
    F_ζ_cor_tracked = [tracked_drifter_data[drifter_num][i].F_ζ_cor for i = 1 : num_iters]
    ζ_visc_tracked = [tracked_drifter_data[drifter_num][i].ζ_visc for i = 1 : num_iters]
    ζ_h_visc_tracked = [tracked_drifter_data[drifter_num][i].ζ_h_visc for i = 1 : num_iters]
    ζ_v_visc_tracked = [tracked_drifter_data[drifter_num][i].ζ_v_visc for i = 1 : num_iters]
    ζ_err_tracked = [tracked_drifter_data[drifter_num][i].ζ_err for i = 1 : num_iters]
    F_ζ_hor_tracked = [tracked_drifter_data[drifter_num][i].F_ζ_hor for i = 1 : num_iters]
    F_ζ_vrt_tracked = [tracked_drifter_data[drifter_num][i].F_ζ_vrt for i = 1 : num_iters]
    ζ_adv_tracked = [tracked_drifter_data[drifter_num][i].ζ_adv for i = 1 : num_iters]
    ζ_h_adv_tracked = [tracked_drifter_data[drifter_num][i].ζ_h_adv for i = 1 : num_iters]
    ζ_interpolated = extract_interpolated_drifter_data(eul_data, "ζ", grid_pos, x, y, t)
    ζ_tendency_interpolated = extract_interpolated_drifter_data(eul_data, "ζ_tendency", grid_pos, x, y, t)
    F_ζ_cor_interpolated = extract_interpolated_drifter_data(eul_data, "F_ζ_cor", grid_pos, x, y, t)
    ζ_visc_interpolated = extract_interpolated_drifter_data(eul_data, "ζ_visc", grid_pos, x, y, t)
    ζ_err_interpolated = extract_interpolated_drifter_data(eul_data, "ζ_err", grid_pos, x, y, t)
    F_ζ_hor_interpolated = extract_interpolated_drifter_data(eul_data, "F_ζ_hor", grid_pos, x, y, t)
    F_ζ_vrt_interpolated = extract_interpolated_drifter_data(eul_data, "F_ζ_vrt", grid_pos, x, y, t)
    ζ_adv_interpolated = extract_interpolated_drifter_data(eul_data, "ζ_adv", grid_pos, x, y, t)
    ζ_h_adv_interpolated = extract_interpolated_drifter_data(eul_data, "ζ_h_adv", grid_pos, x, y, t)

    fig = Figure()
    ax = Axis(fig[1, 1])
    lim = 5e-7
    if plot_mode == "tracked"
        lines!(f*t, ζ_tendency_tracked + ζ_adv_tracked + ζ_err_tracked, label = L"\mathrm{D}\zeta/\mathrm{D}t")
        lines!(f*t, F_ζ_cor_tracked, label = L"\zeta_\text{Cor}")
        lines!(f*t, ζ_v_visc_tracked, label = L"\zeta_\text{visc,v}")
        lines!(f*t, ζ_h_visc_tracked, label = L"\zeta_\text{visc,h}")
        lines!(f*t, F_ζ_hor_tracked, label = L"F_{\zeta,\text{hor}}")
        lines!(f*t, F_ζ_vrt_tracked, label = L"F_{\zeta,\text{vrt}}")
        lines!(f*t, ζ_tendency_tracked + ζ_adv_tracked + ζ_err_tracked - (
            F_ζ_cor_tracked + ζ_v_visc_tracked + ζ_h_visc_tracked + F_ζ_hor_tracked + F_ζ_vrt_tracked),
            label = "residual", color = :black)
        lines!(f*t, ζ_err_tracked, label = L"\zeta_{\text{err}}", color = :black, linestyle = :dot)
        lim = maximum([maximum(abs.(ζ_visc_tracked)), maximum(abs.(F_ζ_hor_tracked))])
    elseif plot_mode == "interpolated"
        lines!(f*t, ζ_tendency_interpolated + ζ_adv_interpolated + ζ_err_interpolated, label = L"\mathrm{D}\zeta/\mathrm{D}t")
        lines!(f*t, F_ζ_cor_interpolated, label = L"\zeta_\text{Cor}")
        lines!(f*t, ζ_visc_interpolated, label = L"\zeta_\text{visc}")
        # lines!(f*t, ζ_err_interpolated, label = L"\zeta_\text{err}")
        lines!(f*t, F_ζ_hor_interpolated, label = L"F_{\zeta,\text{hor}}")
        lines!(f*t, ζ_h_adv_interpolated, label = L"\zeta_\text{adv}")
        # lines!(f*t, F_ζ_vrt_interpolated, label = L"F_{\zeta,\text{vrt}}")
        lines!(f*t, ζ_tendency_interpolated + ζ_adv_interpolated + ζ_err_interpolated - (
            F_ζ_cor_interpolated + ζ_visc_interpolated + F_ζ_hor_interpolated + F_ζ_vrt_interpolated),
            label = "residual")
    end
    ylims!(ax, -lim, lim)
    axislegend(position=:lb)
    display(fig)

end=#

function plot_lagr_ζ_balance(label::String, drifter_num::Int64)

    check_pp_dir(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    i₀ = argmin(abs.(t .- 15/f))
    i₁ = argmin(abs.(t .- 20/f))
    t = t[i₀:i₁]
    tracked_drifter_data = [tracked_drifter_data[n][i₀:i₁] for n in eachindex(tracked_drifter_data)]
    num_iters = length(tracked_drifter_data[drifter_num])

    ζ          = [tracked_drifter_data[drifter_num][i].ζ          for i = 1 : num_iters]
    ζ_tendency = [tracked_drifter_data[drifter_num][i].ζ_tendency for i = 1 : num_iters]
    ζ_adv      = [tracked_drifter_data[drifter_num][i].ζ_adv      for i = 1 : num_iters]
    ζ_err      = [tracked_drifter_data[drifter_num][i].ζ_err      for i = 1 : num_iters]
    F_ζ_hor    = [tracked_drifter_data[drifter_num][i].F_ζ_hor    for i = 1 : num_iters]
    F_ζ_vrt    = [tracked_drifter_data[drifter_num][i].F_ζ_vrt    for i = 1 : num_iters]
    F_ζ_cor    = [tracked_drifter_data[drifter_num][i].F_ζ_cor    for i = 1 : num_iters]
    ζ_h_visc   = [tracked_drifter_data[drifter_num][i].ζ_h_visc   for i = 1 : num_iters]
    ζ_v_visc   = [tracked_drifter_data[drifter_num][i].ζ_v_visc   for i = 1 : num_iters]

    fig = Figure()
    ax = Axis(fig[1, 1])
    lim = 1e-7
    lines!(f*t, ζ_tendency + ζ_adv + ζ_err, label = L"\mathrm{D}\zeta/\mathrm{D}t", color = :black, linestlye = :dash)
    lines!(f*t, F_ζ_hor,  label = L"F_{\zeta,\text{hor}}")
    lines!(f*t, F_ζ_vrt,  label = L"F_{\zeta,\text{vrt}}")
    lines!(f*t, F_ζ_cor,  label = L"\zeta_\text{Cor}")
    lines!(f*t, ζ_h_visc, label = L"\zeta_\text{visc,h}")
    lines!(f*t, ζ_v_visc, label = L"\zeta_\text{visc,v}")
    #=lines!(f*t, ζ_tendency + ζ_adv + ζ_err - (
        F_ζ_cor + ζ_v_visc + ζ_h_visc + F_ζ_hor + F_ζ_vrt),
        label = "residual", color = :black)=#
    lines!(f*t, ζ_err, label = L"\zeta_{\text{err}}", color = :black, linestyle = :dot)
    lines!(f*t, f*ζ, label = L"f\zeta", color = :black, linestyle = :dash)
    lim = maximum([maximum(abs.(ζ_tendency)), maximum(abs.(ζ_adv)),
        maximum(abs.(ζ_err)), maximum(abs.(F_ζ_hor)), maximum(abs.(F_ζ_vrt)),
        maximum(abs.(F_ζ_cor)), maximum(abs.(ζ_h_visc)), maximum(abs.(ζ_v_visc))])
    ylims!(ax, -lim, lim)
    axislegend(position=:lb)
    display(fig)

end

function plot_lagr_δ_balance(label::String, drifter_num::Int64)

    check_pp_dir(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    i₀ = argmin(abs.(t .- 16/f))
    t = t[i₀:end]
    tracked_drifter_data = [tracked_drifter_data[n][i₀:end] for n in eachindex(tracked_drifter_data)]
    num_iters = length(tracked_drifter_data[drifter_num])

    δ          = [tracked_drifter_data[drifter_num][i].δ          for i = 1 : num_iters]
    δ_tendency = [tracked_drifter_data[drifter_num][i].δ_tendency for i = 1 : num_iters]
    δ_adv      = [tracked_drifter_data[drifter_num][i].δ_adv      for i = 1 : num_iters]
    δ_err      = [tracked_drifter_data[drifter_num][i].δ_err      for i = 1 : num_iters]
    F_δ_hor    = [tracked_drifter_data[drifter_num][i].F_δ_hor    for i = 1 : num_iters]
    F_δ_vrt    = [tracked_drifter_data[drifter_num][i].F_δ_vrt    for i = 1 : num_iters]
    F_δ_cor    = [tracked_drifter_data[drifter_num][i].F_δ_cor    for i = 1 : num_iters]
    F_δ_prs    = [tracked_drifter_data[drifter_num][i].F_δ_prs    for i = 1 : num_iters]
    δ_h_visc   = [tracked_drifter_data[drifter_num][i].δ_h_visc   for i = 1 : num_iters]
    δ_v_visc   = [tracked_drifter_data[drifter_num][i].δ_v_visc   for i = 1 : num_iters]

    fig = Figure()
    ax = Axis(fig[1, 1])
    lim = 5e-7
    lines!(f*t, δ_tendency + δ_adv + δ_err, label = L"\mathrm{D}\delta/\mathrm{D}t")
    lines!(f*t, F_δ_hor,  label = L"F_{\delta,\text{hor}}")
    lines!(f*t, F_δ_vrt,  label = L"F_{\delta,\text{vrt}}")
    lines!(f*t, F_δ_cor,  label = L"\delta_\text{Cor}")
    lines!(f*t, δ_h_visc, label = L"\delta_\text{visc,h}")
    lines!(f*t, δ_v_visc, label = L"\delta_\text{visc,v}")
    lines!(f*t, F_δ_prs,  label = L"F_{\delta,\text{prs}}")
    lines!(f*t, δ_tendency + δ_adv + δ_err - (
        F_δ_cor + δ_v_visc + δ_h_visc + F_δ_hor + F_δ_vrt + F_δ_prs),
        label = "residual", color = :black)
    lines!(f*t, δ_err, label = L"\delta_{\text{err}}", color = :black, linestyle = :dot)
    lines!(f*t, f*δ, label = L"f\delta", color = :black, linestyle = :dash)
    lim = maximum([maximum(abs.(δ_tendency)), maximum(abs.(δ_adv)),
        maximum(abs.(δ_err)), maximum(abs.(F_δ_hor)), maximum(abs.(F_δ_vrt)),
        maximum(abs.(F_δ_cor)), maximum(abs.(F_δ_prs)), maximum(abs.(δ_h_visc)),
        maximum(abs.(δ_v_visc))])
    ylims!(ax, -lim, lim)
    axislegend(position=:rb)
    display(fig)

end

investigate_lagr_ζ_balance(label::String) = investigate_lagr_ζ_balance(label, 1)

function investigate_lagr_ζ_balance2(label::String, drifter_num::Int64)
    
    # This compares d/dt(ζ(𝐱(t))) and ζₜ + 𝐮⋅∇ζ

    check_pp_dir(label)
    eul_data = topdata(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    num_iters = length(tracked_drifter_data[drifter_num])

    ζ_tracked = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]
    ζ_tendency_tracked = [tracked_drifter_data[drifter_num][i].ζ_tendency for i = 1 : num_iters]
    ζ_err_tracked = [tracked_drifter_data[drifter_num][i].ζ_err for i = 1 : num_iters]
    ζ_adv_tracked = [tracked_drifter_data[drifter_num][i].ζ_adv for i = 1 : num_iters]

    mean(xs::Vector{Float64}) = (xs[1:end-1] + xs[2:end]) / 2
    Δ(xs::Vector{Float64}) = xs[2:end] - xs[1:end-1]

    Dₜζ_interp = ζ_tendency_tracked + ζ_adv_tracked + ζ_err_tracked
    Δt = Δ(t)
    ∫Dₜζ_interp = [sum([(Dₜζ_interp[j]+Dₜζ_interp[j+1])/2 * Δt[j] for j = 1 : i]) for i = 1 : length(t)-1]
    Dₜζ_lagr = Δ(ζ_tracked) ./ Δt
    ζ̅ = mean(ζ_tracked)
    t̅ = (t[1:end-1] + t[2:end]) / 2

    fig = Figure()
    ax = Axis(fig[1, 1])
    #=lines!(f*t̅, ζ̅)
    lines!(f*t̅, ∫Dₜζ_interp)
    lines!(f*t̅, ∫Dₜζ_interp-ζ̅)
    lim = maximum([maximum(abs.(∫Dₜζ_interp)), maximum(abs.(ζ_tracked))])=#
    #lines!(f*t̅, mean(Dₜζ_interp))
    #lines!(f*t̅, Dₜζ_lagr)
    #lim = maximum([maximum(abs.(Dₜζ_lagr)), maximum(abs.(Dₜζ_interp))])
    lines!(f*t̅, abs.(mean(Dₜζ_interp)) ./ (abs.(mean(Dₜζ_interp)) + abs.(Dₜζ_lagr)))
    lines!(f*t̅, abs.(ζ̅) ./ (abs.(ζ̅) + abs.(∫Dₜζ_interp)))
    lim = 1
    ylims!(ax, -lim, lim)
    display(fig)

end

function isRoughly(x::Float64, y::Float64) :: Bool

    return x * y > 0 && abs(y)/1.5 < abs(x) < 1.5abs(y)

end

function investigate_exceptional_times(drifter_num::Int64)

    eul_data = topdata(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    num_iters = length(tracked_drifter_data[drifter_num])
    iterations = eul_data.iterations

    x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]
    y = [tracked_drifter_data[drifter_num][i].y for i = 1 : num_iters]
    ζ_tracked = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]
    ζ_tendency_tracked = [tracked_drifter_data[drifter_num][i].ζ_tendency for i = 1 : num_iters]
    ζ_cor_tracked = [tracked_drifter_data[drifter_num][i].ζ_cor for i = 1 : num_iters]
    ζ_visc_tracked = [tracked_drifter_data[drifter_num][i].ζ_visc for i = 1 : num_iters]
    ζ_err_tracked = [tracked_drifter_data[drifter_num][i].ζ_err for i = 1 : num_iters]
    F_ζ_hor_tracked = [tracked_drifter_data[drifter_num][i].F_ζ_hor for i = 1 : num_iters]
    F_ζ_vrt_tracked = [tracked_drifter_data[drifter_num][i].F_ζ_vrt for i = 1 : num_iters]
    ζ_adv_tracked = [tracked_drifter_data[drifter_num][i].ζ_adv for i = 1 : num_iters]

    for i = 2 : num_iters - 1
        if !isRoughly(ζ_visc_tracked[i-1], ζ_visc_tracked[i]) &&
                !isRoughly(ζ_visc_tracked[i], ζ_visc_tracked[i+1]) &&
                20 < f*t[i] < 30
            for j = i - 1 : i + 1
                @info j, f*t[j]
            end
            println("")
        end
    end

end
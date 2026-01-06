include("pp-io.jl")
include("drifters-refactored.jl")

label = "just-updated-test"

function f̅(f::Vector{<:AbstractFloat})
    return (f[1:end-1] + f[2:end]) / 2
end

function ddt(f::Vector{<:AbstractFloat}, t::Vector{<:AbstractFloat}; clean = true)
    deriv = (f[1:end-1] - f[2:end]) ./ (t[1:end-1] - t[2:end])
    if ~clean return deriv end
    return [abs.(q) > 1 ? (deriv[i-1] + deriv[i+1])/2 : q for (i, q) in enumerate(deriv)]
end

function isRoughly(x::AbstractFloat, y::AbstractFloat; ε::AbstractFloat = 0.1) :: Bool

    return x * y ≥ 0 && abs(y)/(1+ε) < abs(x) < (1+ε)abs(y)

end

function plot_lagr_ζ_balance(label::String, drifter_num::Int64)

    check_pp_dir(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    #i₀ = argmin(abs.(t .- 15/f))
    #i₁ = argmin(abs.(t .- 20/f))
    i₀ = Int(round(length(tracked_drifter_data[1])/2))
    i₁ = length(tracked_drifter_data[1])
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

    t̅ = f̅(t)
    dζdt = ddt(ζ, t)

    fig = Figure(size=(999,999))
    ax = Axis(fig[1, 1])
    lim = 1e-7
    lines!(f*t, ζ_tendency + ζ_adv, label = L"\mathrm{D}\zeta/\mathrm{D}t", color = :black)
    lines!(f*t̅, dζdt, label = L"\mathrm{D}\zeta/\mathrm{D}t\text{ (Lagr)}")#, color = :black, linestyle = :dash)
    lines!(f*t, F_ζ_hor,  label = L"F_{\zeta,\text{hor}}")
    lines!(f*t, F_ζ_vrt,  label = L"F_{\zeta,\text{vrt}}")
    lines!(f*t, F_ζ_cor,  label = L"\zeta_\text{Cor}")
    lines!(f*t, ζ_h_visc, label = L"\zeta_\text{visc,h}")
    lines!(f*t, ζ_v_visc, label = L"\zeta_\text{visc,v}")
    #=lines!(f*t, ζ_tendency + ζ_adv + ζ_err - (
        F_ζ_cor + ζ_v_visc + ζ_h_visc + F_ζ_hor + F_ζ_vrt),
        label = "residual", color = :black)=#
    lines!(f*t, ζ_err, label = L"\zeta_{\text{err}}", color = :black, linestyle = :dot)
    #lines!(f*t, f*ζ, label = L"f\zeta", color = :black, linestyle = :dash)
    lim = maximum([maximum(abs.(ζ_tendency + ζ_adv + ζ_err)), maximum(abs.(F_ζ_hor)),
                    maximum(abs.(F_ζ_vrt)), maximum(abs.(F_ζ_cor)),
                    maximum(abs.(ζ_h_visc)), maximum(abs.(ζ_v_visc))])
    @info lim
    #ylims!(ax, -lim, lim)
    axislegend(position=:lb)
    display(fig)

end

function plot_lagr_δ_balance(label::String, drifter_num::Int64)

    check_pp_dir(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    i₀ = 300
    i₁ = 500
    t = t[i₀:i₁]
    tracked_drifter_data = [tracked_drifter_data[n][i₀:i₁] for n in eachindex(tracked_drifter_data)]
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

    t̅ = f̅(t)
    dδdt = ddt(δ, t)

    fig = Figure(size=(999,999))
    ax = Axis(fig[1, 1])
    lim = 5e-7
    lines!(f*t, δ_tendency + δ_adv, label = L"\mathrm{D}\delta/\mathrm{D}t")
    lines!(f*t̅, dδdt, label = L"\mathrm{D}\delta/\mathrm{D}t\text{ (Lagr)}")#, color = :black, linestyle = :dash)
    lines!(f*t̅, f̅(δ_tendency + δ_adv) - dδdt)
    lines!(f*t, F_δ_hor,  label = L"F_{\delta,\text{hor}}")
    lines!(f*t, F_δ_vrt,  label = L"F_{\delta,\text{vrt}}")
    lines!(f*t, F_δ_cor,  label = L"\delta_\text{Cor}")
    lines!(f*t, δ_h_visc, label = L"\delta_\text{visc,h}")
    lines!(f*t, δ_v_visc, label = L"\delta_\text{visc,v}")
    lines!(f*t, F_δ_prs,  label = L"F_{\delta,\text{prs}}")
    #=lines!(f*t, δ_tendency + δ_adv + δ_err - (
        F_δ_cor + δ_v_visc + δ_h_visc + F_δ_hor + F_δ_vrt + F_δ_prs),
        label = "residual", color = :black)=#
    lines!(f*t, δ_err, label = L"\delta_{\text{err}}", color = :black, linestyle = :dot)
    lines!(f*t, f*δ, label = L"f\delta", color = :black, linestyle = :dash)
    lim = maximum([maximum(abs.(δ_tendency + δ_adv)),
        maximum(abs.(F_δ_hor)), maximum(abs.(F_δ_vrt)),
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
    tracked_drifter_data = [tracked_drifter_data[i] for i in eachindex(tracked_drifter_data)]
    num_iters = length(tracked_drifter_data[drifter_num])

    x          = [tracked_drifter_data[drifter_num][i].x          for i = 1 : num_iters]
    y          = [tracked_drifter_data[drifter_num][i].y          for i = 1 : num_iters]
    ζ          = [tracked_drifter_data[drifter_num][i].ζ          for i = 1 : num_iters]
    ζ_tendency = [tracked_drifter_data[drifter_num][i].ζ_tendency for i = 1 : num_iters]
    ζ_adv      = [tracked_drifter_data[drifter_num][i].ζ_adv      for i = 1 : num_iters]

    Δ(xs::Vector{<:AbstractFloat}) = xs[2:end] - xs[1:end-1]

    Dₜζ_interp = ζ_tendency + ζ_adv
    Δt = Δ(t)
    ∫Dₜζ_interp = [sum([(Dₜζ_interp[j]+Dₜζ_interp[j+1])/2 * Δt[j] for j = 1 : i]) for i = 1 : length(t)-1]
    Dₜζ_lagr = ddt(ζ, t)
    ζ̅ = f̅(ζ)
    t̅ = f̅(t)

    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(f*t, ζ)
    #=lines!(f*t̅, ζ̅)
    lines!(f*t̅, ∫Dₜζ_interp)
    lines!(f*t̅, ∫Dₜζ_interp-ζ̅)=#
    lim = maximum([maximum(abs.(∫Dₜζ_interp)), maximum(abs.(ζ))])
    #=lines!(f*t̅, mean(Dₜζ_interp))
    lines!(f*t̅, Dₜζ_lagr)
    lim = maximum([maximum(abs.(Dₜζ_lagr)), maximum(abs.(Dₜζ_interp))])=#
    #lines!(f*t̅, abs.(mean(Dₜζ_interp)) ./ (abs.(mean(Dₜζ_interp)) + abs.(Dₜζ_lagr)))
    #lines!(f*t̅, abs.(ζ̅) ./ (abs.(ζ̅) + abs.(∫Dₜζ_interp)))
    #lim = 1
    ylims!(ax, -lim, lim)
    display(fig)

    ζ = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]
    ζ_interp = extract_interpolated_drifter_data(eul_data, "ζ", (Face(), Center()), x, y, t)
    ζ = ζ[1:500]
    ζ_interp = ζ_interp[1:500]
    t = t[1:500]
    t̅ = f̅(t)
    x = x[1:500]
    y = y[1:500]
    fig2 = Figure()
    ax = Axis(fig2[1, 1])
    #scatter!(ax, y[2:end-1], ddt(ddt(ζ, t), t̅))
    #scatter!(ax, y[2:end-1], 3e-3 * ddt(ddt(u_interp, t), t̅))
    lines!(ax, f̅(f̅(t)), ddt(ddt(ζ, t), t̅))
    lines!(ax, f̅(f̅(t)), ddt(ddt(ζ_interp, t), t̅))
    display(fig2)

end

function investigate_exceptional_times(drifter_num::Int64)

    eul_data = topdata(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    num_iters = length(tracked_drifter_data[drifter_num])
    iterations = eul_data.iterations

    x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]
    y = [tracked_drifter_data[drifter_num][i].y for i = 1 : num_iters]
    ζ = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]

    t̅ = f̅(t)
    dζdt = ddt(ζ, t)

    avg_Δt = (t[end] - t[1])/(length(t) - 1)
    Δt = t[2:end] - t[1:end-1]

    for (i, Δ) in enumerate(Δt)
        if !isRoughly(Δ, avg_Δt; ε = 0.1)# && 24 < f*t[i] < 30
            for j = i : i + 1
                @info j, f*t[j]
            end
            println("")
        end
    end

end

function tame_spikes(label, drifter_num)

    check_pp_dir(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    i₀ = 1
    i₁ = length(tracked_drifter_data[1])
    t = t[i₀:i₁]
    tracked_drifter_data = [tracked_drifter_data[n][i₀:i₁] for n in eachindex(tracked_drifter_data)]
    num_iters = length(tracked_drifter_data[drifter_num])

    x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]
    y = [tracked_drifter_data[drifter_num][i].y for i = 1 : num_iters]
    ζ = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]
    ζ_tendency = [tracked_drifter_data[drifter_num][i].ζ_tendency for i = 1 : num_iters]
    ζ_adv      = [tracked_drifter_data[drifter_num][i].ζ_adv      for i = 1 : num_iters]
    Dₜζ_calc = ζ_tendency + ζ_adv

    eul_data = topdata(label)
    ζ_interp = extract_interpolated_drifter_data(eul_data, "ζ", (Face(), Face()), x, y, t)

    t̅ = f̅(t)
    dζdt = ddt(ζ, t)
    dirtydζdt = ddt(ζ, t; clean = false)
    removed_is = filter(1:length(dζdt)) do i dζdt[i] != dirtydζdt[i] end
    dζdt_interp = ddt(ζ_interp, t)

    avg_Δt = (t[end] - t[1])/(length(t) - 1)
    Δt = t[2:end] - t[1:end-1]
    slow_is = filter(1:length(Δt)) do i !isRoughly(Δt[i], avg_Δt; ε = 0.1) end

    Δx = x[2:end] - x[1:end-1]
    xjump_is = filter(1:length(Δx)) do i abs(Δx[i]) > 1e3 end
    Δy = y[2:end] - y[1:end-1]
    yjump_is = filter(1:length(Δy)) do i abs(Δy[i]) > 1e3 end

    fig = Figure(size=(800,500))
    ax = Axis(fig[1, 1])
    lines!(ax, f*t̅, dζdt, label = L"\mathrm{d}\zeta(\mathbf{x}(t))/\mathrm{d}t")
    lines!(ax, f*t̅, dζdt_interp, label = L"\mathrm{d}\zeta(\mathbf{x}(t))/\mathrm{d}t\text{ (PP-tracked)}")
    lines!(ax, f*t, Dₜζ_calc, label = L"\partial_t\zeta+\mathbf{u}\cdot\nabla\zeta")
    scatter!(ax, map(i -> f*t̅[i], slow_is), 0*slow_is, marker = '.', markersize = 30, color = :black)
    scatter!(ax, map(i -> f*t[i], removed_is), 0*removed_is, marker = '+', markersize = 10, color = :red)
    scatter!(ax, map(i -> f*t̅[i], xjump_is), 0*xjump_is, marker = '.', markersize = 30, color = :green)
    scatter!(ax, map(i -> f*t̅[i], yjump_is), 0*yjump_is, marker = '.', markersize = 30, color = :purple)
    axislegend(position=:lb)
    display(fig)

end

function tame_spikes2(label, drifter_num)
    # Now looking at quantities which aren't tracked by drifters except in post-processing

    check_pp_dir(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    i₀ = 1
    i₁ = length(tracked_drifter_data[1])
    t = t[i₀:i₁]
    tracked_drifter_data = [tracked_drifter_data[n][i₀:i₁] for n in eachindex(tracked_drifter_data)]
    num_iters = length(tracked_drifter_data[drifter_num])

    x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]
    y = [tracked_drifter_data[drifter_num][i].y for i = 1 : num_iters]
    ζ = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]

    eul_data = topdata(label)
    u = extract_interpolated_drifter_data(eul_data, "u", (Face(), Center()), x, y, t)
    v = extract_interpolated_drifter_data(eul_data, "v", (Center(), Face()), x, y, t)
    w = extract_interpolated_drifter_data(eul_data, "w", (Center(), Center()), x, y, t)
    b = extract_interpolated_drifter_data(eul_data, "b", (Center(), Center()), x, y, t)
    δ = extract_interpolated_drifter_data(eul_data, "δ", (Center(), Center()), x, y, t)

    t̅ = f̅(t)
    dζdt = ddt(ζ, t)
    dirtydζdt = ddt(ζ, t; clean = false)
    removed_is = filter(1:length(dζdt)) do i dζdt[i] != dirtydζdt[i] end
    dxdt = ddt(x, t; clean = false)
    dydt = ddt(y, t)
    dudt = ddt(u, t)
    dvdt = ddt(v, t)
    dwdt = ddt(w, t)
    dbdt = ddt(b, t)
    dbdt = [abs(q) > 0.00001 ? 0 : q for q in dbdt]
    dδdt = ddt(δ, t)

    avg_Δt = (t[end] - t[1])/(length(t) - 1)
    Δt = t[2:end] - t[1:end-1]
    slow_is = filter(1:length(Δt)) do i !isRoughly(Δt[i], avg_Δt; ε = 0.1) end

    Δx = x[2:end] - x[1:end-1]
    xjump_is = filter(1:length(Δx)) do i abs(Δx[i]) > 1e3 end
    Δy = y[2:end] - y[1:end-1]
    yjump_is = filter(1:length(Δy)) do i abs(Δy[i]) > 1e3 end

    fig = Figure(size=(800,500))
    ax = Axis(fig[1, 1])
    Lx = eul_data.Lx
    lines!(ax, f*t̅, dbdt)
    scatter!(ax, map(i -> f*t̅[i], slow_is), 0*slow_is, marker = '.', markersize = 30, color = :black)
    scatter!(ax, map(i -> f*t[i], removed_is), 0*removed_is, marker = '+', markersize = 10, color = :red)
    scatter!(ax, map(i -> f*t̅[i], xjump_is), 0*xjump_is, marker = '.', markersize = 30, color = :green)
    scatter!(ax, map(i -> f*t̅[i], yjump_is), 0*yjump_is, marker = '.', markersize = 30, color = :purple)
    display(fig)

end

function tame_spikes3()
    label = "just-updated-test"
    drifter_num = 3
    # Now looking at quantities which aren't tracked by drifters except in post-processing

    check_pp_dir(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    i₀ = 1
    i₁ = length(tracked_drifter_data[1])
    t = t[i₀:i₁]
    tracked_drifter_data = [tracked_drifter_data[n][i₀:i₁] for n in eachindex(tracked_drifter_data)]
    num_iters = length(tracked_drifter_data[drifter_num])

    x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]

    eul_data = topdata(label)

    fig = Figure(size=(800,500))
    ax = Axis(fig[1, 1])
    Lx = eul_data.Lx
    lines!(ax, f*t[end-100:end-50], (x[end-100:end-50] .+ 0.5Lx) .% Lx)
    xlims!(ax, 42, 44)
    ylims!(ax, 1e4, 2e4)
    display(fig)

end
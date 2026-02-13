using TickTock
include("pp-io.jl")
include("drifters-refactored.jl")

f = 1e-4
label = "addednhspressure"

indexnearest(arr::Array, val) = argmin(abs.(arr .- val))
indexnearest(range::StepRangeLen, val) = Int(round(Float64((val-range.ref)/range.step))) + range.offset
# Writing this function (rather than using the same as we would for arrays)
# improves things by a factor of roughly 6 with 300 points

function f̅(f::Vector{<:AbstractFloat})
    return (f[1:end-1] + f[2:end]) / 2
end

function ddt(f::Vector{<:AbstractFloat}, t::Vector{<:AbstractFloat}; clean = true)
    deriv = (f[1:end-1] - f[2:end]) ./ (t[1:end-1] - t[2:end])
    if ~clean return deriv end
    return [abs.(q) > 1 ? (deriv[i-1] + deriv[i+1])/2 : q for (i, q) in enumerate(deriv)]
end

Δ(xs::Vector{<:AbstractFloat}) = xs[2:end] - xs[1:end-1]
∫(y_x::Vector{<:AbstractFloat}, x::Vector{<:AbstractFloat}) = [sum([(y_x[j] + y_x[j+1])/2 * Δ(x)[j] for j = 1 : i]) for i = 1 : length(x) - 1]

function isRoughly(x::AbstractFloat, y::AbstractFloat; ε::AbstractFloat = 0.1) :: Bool

    return x * y ≥ 0 && abs(y)/(1+ε) < abs(x) < (1+ε)abs(y)

end

function plot_lagr_ζ_budget(label::String, drifter_num::Int64)

    check_pp_dir(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)

    #i₀ = argmin(abs.(t .- 15/f))
    #i₁ = argmin(abs.(t .- 20/f))
    i₀ = Int(round(length(tracked_drifter_data[1])/2))
    i₁ = length(tracked_drifter_data[1])

    t = t[i₀:i₁]
    drifter = tracked_drifter_data[drifter_num][i₀:i₁]

    ζ        = map(d -> d.ζ,        drifter)
    ζ_t      = map(d -> d.ζ_t,      drifter)
    ζ_adv    = map(d -> d.ζ_adv,    drifter)
    ζ_err    = map(d -> d.ζ_err,    drifter)
    ζ_h_visc = map(d -> d.ζ_h_visc, drifter)
    ζ_v_visc = map(d -> d.ζ_v_visc, drifter)
    F_ζ_cor  = map(d -> d.F_ζ_cor,  drifter)
    F_ζ_hor  = map(d -> d.F_ζ_hor,  drifter)
    # F_ζ_vrt  = map(d -> d.F_ζ_vrt,  drifter)
    # Automtically zero at z = 0 ^

    t̅ = f̅(t)
    dζdt = ddt(ζ, t)

    fig = Figure(size=(999,999))
    ax = Axis(fig[1, 1])
    lim = 1e-7
    lines!(f*t, ζ_t + ζ_adv, label = L"\mathrm{D}\zeta/\mathrm{D}t", color = :black)
    lines!(f*t̅, dζdt, label = L"\mathrm{D}\zeta/\mathrm{D}t\text{ (Lagr)}")#, color = :black, linestyle = :dash)
    lines!(f*t, F_ζ_hor,  label = L"F_{\zeta,\text{hor}}")
    # lines!(f*t, F_ζ_vrt,  label = L"F_{\zeta,\text{vrt}}")
    # Automtically zero ^
    lines!(f*t, F_ζ_cor,  label = L"\zeta_\text{Cor}")
    lines!(f*t, ζ_h_visc, label = L"\zeta_\text{visc,h}")
    lines!(f*t, ζ_v_visc, label = L"\zeta_\text{visc,v}")
    lines!(f*t, ζ_err, label = L"\zeta_{\text{err}}", color = :black, linestyle = :dot)
    lim = maximum([maximum(abs.(ζ_t + ζ_adv + ζ_err)), maximum(abs.(F_ζ_hor)),
                   maximum(abs.(F_ζ_cor)), maximum(abs.(ζ_h_visc)), maximum(abs.(ζ_v_visc))])
    #ylims!(ax, -lim, lim)
    axislegend(position=:lb)
    display(fig)

end

function plot_lagr_δ_budget(label::String, drifter_num::Int64)

    check_pp_dir(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)

    i₀ = argmin(abs.(t .- 5/f))
    i₁ = argmin(abs.(t .- 30/f))
    # i₀ = Int(round(length(tracked_drifter_data[1])/2))
    # i₁ = Int(round(length(tracked_drifter_data[1])*0.8))

    t = t[i₀:i₁]
    drifter = tracked_drifter_data[drifter_num][i₀:i₁]

    δ        = map(d -> d.δ,        drifter)
    δ_t      = map(d -> d.δ_t,      drifter)
    δ_adv    = map(d -> d.δ_adv,    drifter)
    δ_err    = map(d -> d.δ_err,    drifter)
    δ_h_visc = map(d -> d.δ_h_visc, drifter)
    δ_v_visc = map(d -> d.δ_v_visc, drifter)
    F_δ_hor  = map(d -> d.F_δ_hor,  drifter)
    F_δ_vrt  = map(d -> d.F_δ_vrt,  drifter)
    F_δ_cor  = map(d -> d.F_δ_cor,  drifter)
    F_δ_prs  = map(d -> d.F_δ_prs,  drifter)

    t̅ = f̅(t)
    dδdt = ddt(δ, t)

    fig = Figure(size=(999,999))
    ax = Axis(fig[1, 1])
    lim = 1e-7
    lines!(f*t, δ_t + δ_adv, label = L"\mathrm{D}\delta/\mathrm{D}t", color = :black)
    lines!(f*t̅, dδdt, label = L"\mathrm{D}\delta/\mathrm{D}t\text{ (Lagr)}")#, color = :black, linestyle = :dash)
    lines!(f*t̅, f̅(δ_t + δ_adv) - dδdt)
    lines!(f*t, F_δ_hor,  label = L"F_{\delta,\text{hor}}")
    lines!(f*t, F_δ_vrt,  label = L"F_{\delta,\text{vrt}}")
    lines!(f*t, F_δ_cor,  label = L"\delta_\text{Cor}")
    lines!(f*t, δ_h_visc, label = L"\delta_\text{visc,h}")
    lines!(f*t, δ_v_visc, label = L"\delta_\text{visc,v}")
    lines!(f*t, F_δ_prs,  label = L"F_{\delta,\text{prs}}")
    lines!(f*t, δ_t + δ_adv + δ_err - (
        F_δ_cor + δ_v_visc + δ_h_visc + F_δ_hor + F_δ_vrt + F_δ_prs),
        label = "residual", color = :black)
    lines!(f*t, δ_err, label = L"\delta_{\text{err}}", color = :black, linestyle = :dot)
    lines!(f*t, f*δ, label = L"f\delta", color = :black, linestyle = :dash)
    lim = maximum([maximum(abs.(δ_t + δ_adv)),
        maximum(abs.(F_δ_hor)), maximum(abs.(F_δ_vrt)),
        maximum(abs.(F_δ_cor)), maximum(abs.(F_δ_prs)), maximum(abs.(δ_h_visc)),
        maximum(abs.(δ_v_visc))])
    ylims!(ax, -lim, lim)
    axislegend(position=:rb)
    display(fig)

end

function investigate_lagr_balance(varname::String; label::String = "addednhspressure", drifter_num::Int64 = 1)
    # I wrote this function much more intelligently than others
    # This compares d/dt(var(𝐱(t))) and varₜ + 𝐮⋅∇var

    check_pp_dir(label)
    eul_data = topdata(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    # i₀ = argmin(abs.(f*t .- 5))
    # i₁ = argmin(abs.(f*t .- 15))
    i₀ = 1
    i₁ = length(t)

    drifter = tracked_drifter_data[drifter_num][i₀:i₁]
    t = t[i₀:i₁]

    var     = map(d -> d[Symbol(varname)], drifter)
    var_t   = map(d -> d[Symbol(varname * "_t")], drifter)
    var_adv = map(d -> d[Symbol(varname * "_adv")], drifter)

    Dₜvar_eul_interp = var_t + var_adv
    Δt = Δ(t)
    ∫Dₜvar_eul_interp = ∫(Dₜvar_eul_interp, t)
    Dₜvar_lagr = ddt(var, t)
    var̅ = f̅(var)
    t̅ = f̅(t)

    fig = Figure()
    ax = Axis(fig[1, 1])
    offset = var[30] - ∫Dₜvar_eul_interp[30]
    lines!(f*t, var)
    lines!(f*t̅, ∫Dₜvar_eul_interp .+ offset)
    # lines!(f*t̅, f̅(Dₜvar_eul_interp))
    # lines!(f*t̅, Dₜvar_lagr)
    # lines!(f*t, var_t)
    # lines!(f*t, var_adv)
    display(fig)

end
investigate_lagr_ζ_balance() = investigate_lagr_balance("ζ")
investigate_lagr_δ_balance() = investigate_lagr_balance("δ")
# investigate_lagr_u_balance() = investigate_lagr_balance("u")

function highζ_tsample(drifter)

    is = 1 : length(drifter)
    ζ = map(d -> d.ζ, drifter)
    i₀ = argmax(ζ)
    if ζ[i₀] < 10f
        return 1:0, i₀
    end
    i⁻ = argmax(ζ .> f)
    i⁺ = argmax((is .> i₀) .&& (ζ .< f))
    return i⁻:i⁺, i₀

end

function addvarstototaloutput!(varoutput, var, drifter, drifter_bins)
    for (i, d) in enumerate(drifter)
        bin = drifter_bins[i]
        bin = bin == 0 ? 1 : bin
        varoutput[bin] += var.func(d)
    end
end

function budgetaggregated(vars, label = "more-drifters", n_bins = 300)

    t, drifters = extract_tracked_drifter_data(label)
    n_drifters = length(drifters)
    n_vars = length(vars)

    ranges = UnitRange{Int64}[]
    centres = Int64[]
    minrelativet = 0.0
    maxrelativet = 0.0
    for drifter in drifters
        r, c = highζ_tsample(drifter)
        push!(ranges, r)
        push!(centres, c)
        if length(r) > 0
            minrelativet = minimum([minrelativet, t[r[1]]-t[c]])
            maxrelativet = maximum([maxrelativet, t[r[end]]-t[c]])
        end
    end
    Δt = (maxrelativet - minrelativet) / n_bins
    relativetbins = range(start = minrelativet + Δt/2, stop = maxrelativet - Δt/2, length = n_bins)
    varoutputstotal = [zeros(n_bins) for i = 1 : n_vars]
    freq = zeros(n_bins)

    tick()
    for i = 1 : n_drifters
        relativets = t[ranges[i]] .- t[centres[i]]
        drifter_bins = [indexnearest(relativetbins, relt) for relt in relativets]
        for bin in drifter_bins
            bin = bin == 0 ? 1 : bin
            freq[bin] += 1
        end
        for j = 1 : n_vars
            addvarstototaloutput!(varoutputstotal[j], vars[j], drifters[i][ranges[i]], drifter_bins)
        end
    end
    @info tok()
    varoutputsmean = [varoutputtotal ./ freq for varoutputtotal in varoutputstotal]

    fig = Figure()
    ax = Axis(fig[1, 1])
    for (i, var) in enumerate(vars)
        lines!(ax, f*relativetbins, varoutputsmean[i], label = var.label)
    end
    axislegend(position=:lb)
    display(fig)

end

ζbudgetvars  = ((func = d -> d.ζ_t + d.ζ_adv,       label = L"\mathrm{D}\zeta/\mathrm{D}t"),
                (func = d -> d.F_ζ_cor,             label = L"F_{\zeta\text{, Cor}}"),
                (func = d -> d.F_ζ_hor,             label = L"F_{\zeta\text{, hor}}"),
                (func = d -> d.ζ_h_visc,            label = L"\zeta_\text{visc, h}"),
                (func = d -> d.ζ_v_visc,            label = L"\zeta_\text{visc, v}"))
δbudgetvars  = ((func = d -> d.δ_t + d.δ_adv,       label = L"\mathrm{D}\delta/\mathrm{D}t"),
                (func = d -> d.F_δ_cor,             label = L"F_{\delta\text{, Cor}}"),
                (func = d -> d.F_δ_hor,             label = L"F_{\delta\text{, hor}}"),
                (func = d -> d.F_δ_prs,             label = L"F_{\delta\text{, prs}}"),
                (func = d -> d.δ_h_visc,            label = L"\delta_\text{visc, h}"),
                (func = d -> d.δ_v_visc,            label = L"\delta_\text{visc, v}"))
δbudgetvars2 = ((func = d -> d.δ_t + d.δ_adv,       label = L"\mathrm{D}\delta/\mathrm{D}t"),
                (func = d -> d.F_δ_cor + d.F_δ_prs, label = L"F_{\delta\text{, Cor}}+F_{\delta\text{, prs}}"),
                (func = d -> d.F_δ_hor,             label = L"F_{\delta\text{, hor}}"),
                (func = d -> d.δ_h_visc,            label = L"\delta_\text{visc, h}"),
                (func = d -> d.δ_v_visc,            label = L"\delta_\text{visc, v}"))
testbudg() = budgetaggregated(testplotvars)
ζbudg() = budgetaggregated(ζbudgetvars)
δbudg() = budgetaggregated(δbudgetvars)
δbudg2() = budgetaggregated(δbudgetvars2)

ζbudg()

# function test_eul_var(varname::String)
# 
#     # This tests whether expected and actualy derivatives match
# 
#     label = "addednhspressure"
#     data = topdata(label)
#     iterations = data.iterations
#     file = data.file
#     var        = [file["timeseries/$varname/$iter"][1, 1, 1]          for iter in iterations]
#     t          = [file["timeseries/t/$iter"]                          for iter in iterations]
#     var_t_pred = [file["timeseries/" * varname * "_t/$iter"][1, 1, 1] for iter in iterations]
#     Δ(vec) = vec[1:end-1] - vec[2:end]
#     f̅(vec) = (vec[1:end-1] + vec[2:end]) / 2
#     var_t_actual = Δ(var) ./ Δ(t)
#     t̅ = f̅(t)
# 
#     fig = Figure()
#     ax = Axis(fig[1, 1])
#     lines!(ax, t, var_t_pred)
#     lines!(ax, t̅, var_t_actual)
#     display(fig)
# 
# end
# test_eul_δ() = test_eul_var("δ")
# test_eul_u() = test_eul_var("u")
# test_eul_ζ() = test_eul_var("ζ")
#
################################################
## Some functions comparing particle tracking ##
##   and post-processing interpolation onto   ##
##             particle positions             ##
################################################
#
# function test_δ_interpolation()
# 
#     f = 1e-4
#     label = "addednhspressure"
#     t_full, tracked_drifter_data = extract_tracked_drifter_data(label)
#     ft_full = f*t_full
# 
#     i₀ = argmin(abs.(ft_full .- 5))
#     i₁ = argmin(abs.(ft_full .- 15))
#     drifter_full = tracked_drifter_data[1]
#     num_iters_full = length(drifter_full)
#     drifter = drifter_full[i₀:i₁]
#     t = t_full[i₀:i₁]
#     ft = ft_full[i₀:i₁]
#     num_iters = length(drifter)
#     x_full      = [drifter_full[i].x      for i = 1 : num_iters_full]
#     y_full      = [drifter_full[i].y      for i = 1 : num_iters_full]
#     x           = [drifter[i].x           for i = 1 : num_iters]
#     y           = [drifter[i].y           for i = 1 : num_iters]
#     δ           = [drifter[i].δ           for i = 1 : num_iters]
#     δ_t  = [drifter[i].δ_t  for i = 1 : num_iters]
#     δ_adv       = [drifter[i].δ_adv       for i = 1 : num_iters]
#     eul_data = topdata(label)
#     @info length(t)
#     @info length(t_full)
#     δ2 = extract_interpolated_drifter_data(eul_data, "δ", (Center(), Center()), x_full, y_full, t_full)
#     δ_t2 = extract_interpolated_drifter_data(eul_data, "δ_t", (Center(), Center()), x_full, y_full, t_full)
#     δ_adv2 = extract_interpolated_drifter_data(eul_data, "δ_adv", (Center(), Center()), x_full, y_full, t_full)
#     δ2 = δ2[i₀:i₁]
#     δ_t2 = δ_t2[i₀:i₁]
#     δ_adv2 = δ_adv2[i₀:i₁]
# 
#     fig = Figure()
#     ax = Axis(fig[1, 1])
#     lines!(ax, ft, δ, label = L"\delta")
#     lines!(ax, ft, δ2, label = L"\delta_I")
#     lines!(ax, ft, δ_t/f, label = L"\delta_t/f")
#     lines!(ax, ft, δ_t2/f, label = L"\delta_{t,I}/f")
#     lines!(ax, ft, δ_adv/f, label = L"\mathbf{u}\cdot\nabla\delta/f")
#     lines!(ax, ft, δ_adv2/f, label = L"[\mathbf{u}\cdot\nabla\delta]_I/f")
#     axislegend(position=:lb)
#     display(fig)
# 
# end
#
# function test_u_interpolation()
# 
#     f = 1e-4
#     label = "addednhspressure"
#     t_full, tracked_drifter_data = extract_tracked_drifter_data(label)
#     ft_full = f*t_full
# 
#     # i₀ = argmin(abs.(ft_full .- 5))
#     # i₁ = argmin(abs.(ft_full .- 25))
#     i₀ = 1
#     i₁ = length(tracked_drifter_data[1])
#     drifter_full = tracked_drifter_data[1]
#     num_iters_full = length(drifter_full)
#     drifter = drifter_full[i₀:i₁]
#     t = t_full[i₀:i₁]
#     ft = ft_full[i₀:i₁]
#     num_iters = length(drifter)
#     x_full      = [drifter_full[i].x      for i = 1 : num_iters_full]
#     y_full      = [drifter_full[i].y      for i = 1 : num_iters_full]
#     x           = [drifter[i].x           for i = 1 : num_iters]
#     y           = [drifter[i].y           for i = 1 : num_iters]
#     u           = [drifter[i].u           for i = 1 : num_iters] .- 0.5
#     u_t  = [drifter[i].u_t  for i = 1 : num_iters]
#     u_adv       = [drifter[i].u_adv       for i = 1 : num_iters]
#     eul_data = topdata(label)
#     @info length(t)
#     @info length(t_full)
#     u2 = extract_interpolated_drifter_data(eul_data, "u", (Center(), Center()), x_full, y_full, t_full)
#     u_t2 = extract_interpolated_drifter_data(eul_data, "u_t", (Center(), Center()), x_full, y_full, t_full)
#     u_adv2 = extract_interpolated_drifter_data(eul_data, "u_adv", (Center(), Center()), x_full, y_full, t_full)
#     u2 = u2[i₀:i₁]
#     u_t2 = u_t2[i₀:i₁]
#     u_adv2 = u_adv2[i₀:i₁]
# 
#     fig = Figure()
#     ax = Axis(fig[1, 1])
#     lines!(ax, ft, u, label = L"u")
#     lines!(ax, ft, u2, label = L"u_I")
#     lines!(ax, ft, u_t/f, label = L"u_t/f")
#     lines!(ax, ft, u_t2/f, label = L"u_{t,I}/f")
#     lines!(ax, ft, u_adv/f, label = L"\mathbf{u}\cdot\nabla u/f")
#     lines!(ax, ft, u_adv2/f, label = L"[\mathbf{u}\cdot\nabla u]_I/f")
#     axislegend(position=:lb)
#     display(fig)
# 
# end
#
# function test_ζ_interpolation()
# 
#     f = 1e-4
#     label = "addednhspressure"
#     t_full, tracked_drifter_data = extract_tracked_drifter_data(label)
#     ft_full = f*t_full
# 
#     i₀ = argmin(abs.(ft_full .- 5))
#     i₁ = argmin(abs.(ft_full .- 15))
#     drifter_full = tracked_drifter_data[1]
#     num_iters_full = length(drifter_full)
#     drifter = drifter_full[i₀:i₁]
#     t = t_full[i₀:i₁]
#     ft = ft_full[i₀:i₁]
#     num_iters = length(drifter)
#     x_full      = [drifter_full[i].x      for i = 1 : num_iters_full]
#     y_full      = [drifter_full[i].y      for i = 1 : num_iters_full]
#     x           = [drifter[i].x           for i = 1 : num_iters]
#     y           = [drifter[i].y           for i = 1 : num_iters]
#     ζ           = [drifter[i].ζ           for i = 1 : num_iters]
#     eul_data = topdata(label)
#     @info length(t)
#     @info length(t_full)
#     ζ2 = extract_interpolated_drifter_data(eul_data, "ζ", (Center(), Center()), x_full, y_full, t_full)
#     ζ2 = ζ2[i₀:i₁]
# 
#     fig = Figure()
#     ax = Axis(fig[1, 1])
#     lines!(ax, ft, ζ, label = L"\zeta")
#     lines!(ax, ft, ζ2, label = L"\zeta_I")
#     axislegend(position=:lb)
#     display(fig)
# 
# end
#
# BELOW ARE SOME FUNCTIONS FOR LOOKING AT WHERE DATA WAS SPIKY IN THE PAST
# I THINK THAT ALL OF THIS HAS NOW BEEN FIXED
#
# function investigate_exceptional_times(drifter_num::Int64)
# 
#     eul_data = topdata(label)
#     t, tracked_drifter_data = extract_tracked_drifter_data(label)
#     num_iters = length(tracked_drifter_data[drifter_num])
#     iterations = eul_data.iterations
# 
#     x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]
#     y = [tracked_drifter_data[drifter_num][i].y for i = 1 : num_iters]
#     ζ = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]
# 
#     t̅ = f̅(t)
#     dζdt = ddt(ζ, t)
# 
#     avg_Δt = (t[end] - t[1])/(length(t) - 1)
#     Δt = t[2:end] - t[1:end-1]
# 
#     for (i, Δ) in enumerate(Δt)
#         if !isRoughly(Δ, avg_Δt; ε = 0.1)# && 24 < f*t[i] < 30
#             for j = i : i + 1
#                 @info j, f*t[j]
#             end
#             println("")
#         end
#     end
# 
# end
# 
# function tame_spikes(label, drifter_num)
# 
#     check_pp_dir(label)
#     t, tracked_drifter_data = extract_tracked_drifter_data(label)
#     i₀ = 1
#     i₁ = length(tracked_drifter_data[1])
#     t = t[i₀:i₁]
#     tracked_drifter_data = [tracked_drifter_data[n][i₀:i₁] for n in eachindex(tracked_drifter_data)]
#     num_iters = length(tracked_drifter_data[drifter_num])
# 
#     x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]
#     y = [tracked_drifter_data[drifter_num][i].y for i = 1 : num_iters]
#     ζ = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]
#     ζ_t   = [tracked_drifter_data[drifter_num][i].ζ_t   for i = 1 : num_iters]
#     ζ_adv = [tracked_drifter_data[drifter_num][i].ζ_adv for i = 1 : num_iters]
#     Dₜζ_calc = ζ_t + ζ_adv
# 
#     eul_data = topdata(label)
#     ζ_interp = extract_interpolated_drifter_data(eul_data, "ζ", (Face(), Face()), x, y, t)
# 
#     t̅ = f̅(t)
#     dζdt = ddt(ζ, t)
#     dirtydζdt = ddt(ζ, t; clean = false)
#     removed_is = filter(1:length(dζdt)) do i dζdt[i] != dirtydζdt[i] end
#     dζdt_interp = ddt(ζ_interp, t)
# 
#     avg_Δt = (t[end] - t[1])/(length(t) - 1)
#     Δt = t[2:end] - t[1:end-1]
#     slow_is = filter(1:length(Δt)) do i !isRoughly(Δt[i], avg_Δt; ε = 0.1) end
# 
#     Δx = x[2:end] - x[1:end-1]
#     xjump_is = filter(1:length(Δx)) do i abs(Δx[i]) > 1e3 end
#     Δy = y[2:end] - y[1:end-1]
#     yjump_is = filter(1:length(Δy)) do i abs(Δy[i]) > 1e3 end
# 
#     fig = Figure(size=(800,500))
#     ax = Axis(fig[1, 1])
#     lines!(ax, f*t̅, dζdt, label = L"\mathrm{d}\zeta(\mathbf{x}(t))/\mathrm{d}t")
#     lines!(ax, f*t̅, dζdt_interp, label = L"\mathrm{d}\zeta(\mathbf{x}(t))/\mathrm{d}t\text{ (PP-tracked)}")
#     lines!(ax, f*t, Dₜζ_calc, label = L"\partial_t\zeta+\mathbf{u}\cdot\nabla\zeta")
#     scatter!(ax, map(i -> f*t̅[i], slow_is), 0*slow_is, marker = '.', markersize = 30, color = :black)
#     scatter!(ax, map(i -> f*t[i], removed_is), 0*removed_is, marker = '+', markersize = 10, color = :red)
#     scatter!(ax, map(i -> f*t̅[i], xjump_is), 0*xjump_is, marker = '.', markersize = 30, color = :green)
#     scatter!(ax, map(i -> f*t̅[i], yjump_is), 0*yjump_is, marker = '.', markersize = 30, color = :purple)
#     axislegend(position=:lb)
#     display(fig)
# 
# end
# 
# function tame_spikes2(label, drifter_num)
#     # Now looking at quantities which aren't tracked by drifters except in post-processing
# 
#     check_pp_dir(label)
#     t, tracked_drifter_data = extract_tracked_drifter_data(label)
#     i₀ = 1
#     i₁ = length(tracked_drifter_data[1])
#     t = t[i₀:i₁]
#     tracked_drifter_data = [tracked_drifter_data[n][i₀:i₁] for n in eachindex(tracked_drifter_data)]
#     num_iters = length(tracked_drifter_data[drifter_num])
# 
#     x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]
#     y = [tracked_drifter_data[drifter_num][i].y for i = 1 : num_iters]
#     ζ = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]
# 
#     eul_data = topdata(label)
#     u = extract_interpolated_drifter_data(eul_data, "u", (Face(), Center()), x, y, t)
#     v = extract_interpolated_drifter_data(eul_data, "v", (Center(), Face()), x, y, t)
#     w = extract_interpolated_drifter_data(eul_data, "w", (Center(), Center()), x, y, t)
#     b = extract_interpolated_drifter_data(eul_data, "b", (Center(), Center()), x, y, t)
#     δ = extract_interpolated_drifter_data(eul_data, "δ", (Center(), Center()), x, y, t)
# 
#     t̅ = f̅(t)
#     dζdt = ddt(ζ, t)
#     dirtydζdt = ddt(ζ, t; clean = false)
#     removed_is = filter(1:length(dζdt)) do i dζdt[i] != dirtydζdt[i] end
#     dxdt = ddt(x, t; clean = false)
#     dydt = ddt(y, t)
#     dudt = ddt(u, t)
#     dvdt = ddt(v, t)
#     dwdt = ddt(w, t)
#     dbdt = ddt(b, t)
#     dbdt = [abs(q) > 0.00001 ? 0 : q for q in dbdt]
#     dδdt = ddt(δ, t)
# 
#     avg_Δt = (t[end] - t[1])/(length(t) - 1)
#     Δt = t[2:end] - t[1:end-1]
#     slow_is = filter(1:length(Δt)) do i !isRoughly(Δt[i], avg_Δt; ε = 0.1) end
# 
#     Δx = x[2:end] - x[1:end-1]
#     xjump_is = filter(1:length(Δx)) do i abs(Δx[i]) > 1e3 end
#     Δy = y[2:end] - y[1:end-1]
#     yjump_is = filter(1:length(Δy)) do i abs(Δy[i]) > 1e3 end
# 
#     fig = Figure(size=(800,500))
#     ax = Axis(fig[1, 1])
#     Lx = eul_data.Lx
#     lines!(ax, f*t̅, dbdt)
#     scatter!(ax, map(i -> f*t̅[i], slow_is), 0*slow_is, marker = '.', markersize = 30, color = :black)
#     scatter!(ax, map(i -> f*t[i], removed_is), 0*removed_is, marker = '+', markersize = 10, color = :red)
#     scatter!(ax, map(i -> f*t̅[i], xjump_is), 0*xjump_is, marker = '.', markersize = 30, color = :green)
#     scatter!(ax, map(i -> f*t̅[i], yjump_is), 0*yjump_is, marker = '.', markersize = 30, color = :purple)
#     display(fig)
# 
# end
# 
# function tame_spikes3()
#     label = "whatisgoingonwithdelta"
#     drifter_num = 3
#     # Now looking at quantities which aren't tracked by drifters except in post-processing
# 
#     check_pp_dir(label)
#     t, tracked_drifter_data = extract_tracked_drifter_data(label)
#     i₀ = 1
#     i₁ = length(tracked_drifter_data[1])
#     t = t[i₀:i₁]
#     tracked_drifter_data = [tracked_drifter_data[n][i₀:i₁] for n in eachindex(tracked_drifter_data)]
#     num_iters = length(tracked_drifter_data[drifter_num])
# 
#     x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]
# 
#     eul_data = topdata(label)
# 
#     fig = Figure(size=(800,500))
#     ax = Axis(fig[1, 1])
#     Lx = eul_data.Lx
#     lines!(ax, f*t[end-100:end-50], (x[end-100:end-50] .+ 0.5Lx) .% Lx)
#     xlims!(ax, 42, 44)
#     ylims!(ax, 1e4, 2e4)
#     display(fig)
# 
# end
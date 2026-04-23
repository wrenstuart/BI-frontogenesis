using JLD2
using CairoMakie
include("pp-io.jl")
include("drifters-refactored.jl")

ζ(drifter) = map(d -> d.ζ, drifter)
δ(drifter) = map(d -> d.δ, drifter)

f = 1e-4

function addphasepath!(t_reptot, n_visitstot, xs, ys, t::Vector{<:AbstractFloat}, drifter, var1::Symbol, var2::Symbol, t_repfunc::Function)
    n_iters = length(drifter)
    for i = 1 : n_iters
        i₁ = argmin(abs.(xs .- drifter[i][var1]))
        i₂ = argmin(abs.(ys .- drifter[i][var2]))
        # I can probably get this to work faster using rounding, given that I
        # know that xs and ys will be nice ranges
        t_reptot[i₁, i₂] += t_repfunc(t[i])
        n_visitstot[i₁, i₂] += 1
    end
end

function t_repfuncgenerate_ζmax(drifter, t)
    # returns a function that maps the time to a value representing how close
    # this drifter's path is to being at maximal ζ
    tζmax = t[argmax(ζ(drifter))]
    t₀ = t[1]
    t₁ = t[end]
    t_span = max(tζmax-t₀, t₁-tζmax)
    if t_span == 0
        @info "aaa"
    end
    return t -> (t-tζmax)/t_span
end

function phasepaths(label::String, var1::Symbol = :ζ, var2::Symbol = :δ, t_repfuncgenerate::Function = t_repfuncgenerate_ζmax, n_bins::Int = 400)

    t, drifters = extract_tracked_drifter_data(label)

    max_x = maximum(map(drifter -> maximum(map(d -> d[var1], drifter)), drifters))
    max_y = maximum(map(drifter -> maximum(map(d -> d[var2], drifter)), drifters))
    min_x = minimum(map(drifter -> minimum(map(d -> d[var1], drifter)), drifters))
    min_y = minimum(map(drifter -> minimum(map(d -> d[var2], drifter)), drifters))
    Δx = (max_x - min_x) / n_bins
    Δy = (max_y - min_y) / n_bins
    xs = range(start = min_x + Δx/2, stop = max_x - Δx/2, length = n_bins)
    ys = range(start = min_y + Δy/2, stop = max_y - Δy/2, length = n_bins)

    t_reptot = zeros(n_bins, n_bins)
    n_visitstot = zeros(n_bins, n_bins)
    for drifter in drifters
        t_repfunc = t_repfuncgenerate(drifter, t)
        addphasepath!(t_reptot, n_visitstot, xs, ys, t, drifter, var1, var2, t_repfunc)
    end
    t_rep = t_reptot ./ n_visitstot
    # NaN renders as transparent, which is what we want

    fig = Figure()
    ax = Axis(fig[1, 1][1, 1], xlabel = L"$\zeta/f$", ylabel = L"$\delta/f$")
    hm = heatmap!(ax, xs/f, ys/f, t_rep; colormap = :seismic, colorrange = (-1, 1))
    Colorbar(fig[1, 1][1, 2], hm)
    display(fig)

end
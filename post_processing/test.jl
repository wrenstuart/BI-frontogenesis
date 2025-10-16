using CairoMakie
using Oceananigans
using JLD2
using Printf
using Oceananigans.Units
using Makie.Colors
using Statistics
using OffsetArrays
using FFTW

include("drifters-refactored.jl")
include("pp-io.jl")

# Set the two dimensional parameters
H = 50    # Depth of mixed layer
f = 1e-4  # Coriolis parameter

function ft_display(t::Float64)
    ft = @sprintf("%.2f", f*t)
    return L"ft=%$(ft)"
end

function ani_xy(label::String, a::Float64, b::Float64)  # Animate b, ζ and δ at surface

    check_pp_dir(label)
    filename_xy_top = data_dir(label) * "BI_xy.jld2"

    # Read in the first iteration. We do this to load the grid
    b_ic = FieldTimeSeries(filename_xy_top, "b", iterations = 0)
    ζ_ic = FieldTimeSeries(filename_xy_top, "ζ", iterations = 0)
    δ_ic = FieldTimeSeries(filename_xy_top, "δ", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xb, yb, ~ = nodes(b_ic)
    xζ, yζ, ~ = nodes(ζ_ic)
    xδ, yδ, ~ = nodes(δ_ic)

    # Now, open the file with our data
    file = jldopen(filename_xy_top)

    # Extract the values that iter can take
    iterations = parse.(Int, keys(file["timeseries/t"]))

    # Set up observables for plotting that will update as the iteration number is updated
    iter = 0
    iter_index = findfirst(x -> x == iter, iterations)
    ζ_xy = file["timeseries/ζ/$iter"][:, :, 1]
    ζ_on_f = ζ_xy/f
    # δ = lift(iter -> file["timeseries/δ/$iter"][:, :, 1], iter)
    δ_on_f = lift(iter -> file["timeseries/δ/$iter"][:, :, 1]/f, iter)
    b_xy = lift(iter -> file["timeseries/b/$iter"][:, :, 1], iter)
    t = lift(iter -> file["timeseries/t/$iter"], iter)   # Time elapsed by this iteration
    str_ft = lift(t -> ft_display(t), t)

    # Calculate the maximum relative vorticity and buoyancy to set the scale for the colourmap
    b_max = maximum(b_ic)
    ζ_max = 20f

    ~, drifters = extract_tracked_drifter_data(label)
    @info length(t)
    @info length(iterations)
    xs = [drifters[i][1].x for i in eachindex(drifters)]
    ys = [drifters[i][1].y for i in eachindex(drifters)]

    @info "Drawing first frame"

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    fig = Figure(size = (950, 320))
    ax_b = Axis(fig[1, 1][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Buoyancy, }b", width = 200, height = 200)
    ax_ζ = Axis(fig[1, 2][1, 1], xlabel = L"$x/\mathrm{km}$", title = L"\text{Vertical vorticity, }\zeta/f", width = 200, height = 200)
    ax_δ = Axis(fig[1, 3][1, 1], xlabel = L"$x/\mathrm{km}$", title = L"\text{Horizontal divergence, }\delta/f", width = 200, height = 200)
    hm_b = heatmap!(ax_b, xb/1kilometer, yb/1kilometer, b_xy; colorrange = (-0.5*b_max, 1.5*b_max));
    sc_b = scatter!(ax_b, xs, ys, marker = '.', markersize = 15, color = :black, transparency = true)
    hm_ζ = heatmap!(ax_ζ, xζ/1kilometer, yζ/1kilometer, ζ_on_f; colormap = :coolwarm, colorrange = (-ζ_max/f, ζ_max/f))
    hm_δ = heatmap!(ax_δ, xδ/1kilometer, yδ/1kilometer, δ_on_f; colormap = :coolwarm, colorrange = (-ζ_max/f, ζ_max/f))
    Colorbar(fig[1, 1][1, 2], hm_b, height = 200)
    Colorbar(fig[1, 2][1, 2], hm_ζ, height = 200)
    Colorbar(fig[1, 3][1, 2], hm_δ, height = 200)
    Makie.Label(fig[0, 1:3], str_ft)

    display(fig)
    
    #=@info "Making an animation from saved data..."
    CairoMakie.record(i -> iter[] = i,
           fig,
           pp_dir(label) * "bζδ-top-vid.mp4",
           iterations[maximum([Int64(round(length(iterations) * a )), 1]) : 2 : Int64(round(length(iterations) * b))],
           framerate = 20)=#

end

# ani_xy("big_test", 0.0, 0.3)

label = "big_test"

check_pp_dir(label)
filename_xy_top = data_dir(label) * "BI_xy.jld2"

# Read in the first iteration. We do this to load the grid
b_ic = FieldTimeSeries(filename_xy_top, "b", iterations = 0)
ζ_ic = FieldTimeSeries(filename_xy_top, "ζ", iterations = 0)
δ_ic = FieldTimeSeries(filename_xy_top, "δ", iterations = 0)

# Load in co-ordinate arrays
# We do this separately for each variable since Oceananigans uses a staggered grid
xb, yb, ~ = nodes(b_ic)

fig = Figure(size = (950, 320))
ax = Axis(fig[1, 1][1, 1], width = 200, height = 200)
x = 1.0 : 10.0
y = 1.0 : 10.0
hm = heatmap!(ax, xb/1kilometer, yb/1kilometer, (x, y) -> sin(x+y))
scatter!(ax, randn(10) .+ 5, randn(10) .+ 5, marker = '.', markersize = 50, color = :red)
Colorbar(fig[1, 1][1, 2], hm, height = 200)
display(fig)
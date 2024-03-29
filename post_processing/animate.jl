using CairoMakie
using Oceananigans
using JLD2
using Printf
using Oceananigans.Units
using Makie.Colors

function ani_xy(label)

    # Set the two dimensional parameters
    H = 50    # Depth of mixed layer
    f = 1e-4  # Coriolis parameter

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"

    # Read in the first iteration. We do this to load the grid
    b_ic = FieldTimeSeries(filename_xy_top, "b", iterations = 0)
    ζ₃_ic = FieldTimeSeries(filename_xy_top, "ζ₃", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xb, yb, zb = nodes(b_ic)
    xζ₃, yζ₃, zζ₃ = nodes(ζ₃_ic)

    # Now, open the file with our data
    file = jldopen(filename_xy_top)

    # Set up observables for plotting that will update as the iteration number is updated
    iter = Observable(0)
    ζ₃_xy = lift(iter -> file["timeseries/ζ₃/$iter"][:, :, 1], iter) # Surface vertical vorticity at this iteration
    ζ_on_f = lift(iter -> ζ₃_xy[]/f, iter)
    b_xy = lift(iter -> file["timeseries/b/$iter"][:, :, 1], iter)   # Surface buoyancy at this iteration
    t = lift(iter -> file["timeseries/t/$iter"], iter)   # Time elapsed by this iteration

    # Extract the values that iter can take
    iterations = parse.(Int, keys(file["timeseries/t"]))

    #=
    t_save = zeros(length(iterations))
    # This will contain the actual physical time elapsed
    =#

    # Calculate the maximum relative vorticity and buoyancy flux to set the scale for the colourmap
    ζ₃_max = 0
    b_max = maximum(b_ic)

    for i = 11 : length(iterations)
        iter[] = iterations[i]
        ζ₃_max = maximum([ζ₃_max, maximum(ζ₃_xy[])])
    end

    ζ₃_max = minimum([ζ₃_max, 20*f])

    @info "Drawing first frame"

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    iter[] = 0
    fig = Figure(size = (1280, 720))
    ax_b = Axis(fig[1, 1][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Buoyancy, }b")
    ax_ζ = Axis(fig[1, 2][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Vertical vorticity, }\zeta/f")
    hm_b = heatmap!(ax_b, xb/1kilometer, yb/1kilometer, b_xy; colorrange = (-0.5*b_max, 1.5*b_max));
    hm_ζ₃ = heatmap!(ax_ζ, xζ₃/1kilometer, yζ₃/1kilometer, ζ_on_f; colormap = :coolwarm, colorrange = (-ζ₃_max/f, ζ₃_max/f));
    Colorbar(fig[1, 1][1, 2], hm_b)
    Colorbar(fig[1, 2][1, 2], hm_ζ₃)

    ##############################################
    # SEE IF CAN CHANGE COLORS, ADD TIME COUNTER #
    ##############################################

    display(fig)
    
    @info "Making an animation from saved data..."
    record(i -> iter[] = i,
           fig,
           "pretty_things/" * label * ".mp4",
           iterations[Int64(round(length(iterations)*0.5)) : length(iterations)],
           framerate = 20)

end

function ani_xz(label)

    # Set the two dimensional parameters
    H = 50    # Depth of mixed layer
    f = 1e-4  # Coriolis parameter

    filename_xz_top = "raw_data/" * label * "_BI_xz" * ".jld2"

    # Read in the first iteration. We do this to load the grid
    b_ic = FieldTimeSeries(filename_xz_top, "b", iterations = 0)
    ζ₃_ic = FieldTimeSeries(filename_xz_top, "ζ₃", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xb, yb, zb = nodes(b_ic)
    xζ₃, yζ₃, zζ₃ = nodes(ζ₃_ic)

    # Now, open the file with our data
    file = jldopen(filename_xz_top)

    # Set up observables for plotting that will update as the iteration number is updated
    iter = Observable(0)
    ζ₃_xz = lift(iter -> file["timeseries/ζ₃/$iter"][:, 1, :], iter) # Surface vertical vorticity at this iteration
    ζ_on_f = lift(iter -> ζ₃_xz[]/f, iter)
    b_xz = lift(iter -> file["timeseries/b/$iter"][:, 1, :], iter)   # Surface buoyancy at this iteration
    t = lift(iter -> file["timeseries/t/$iter"], iter)   # Time elapsed by this iteration

    # Extract the values that iter can take
    iterations = parse.(Int, keys(file["timeseries/t"]))

    #=
    t_save = zeros(length(iterations))
    # This will contain the actual physical time elapsed
    =#

    # Calculate the maximum relative vorticity and buoyancy flux to set the scale for the colourmap
    ζ₃_max = 0
    b_max = maximum([maximum(b_ic), maximum(-b_ic)])

    for i = 11 : length(iterations)
        iter[] = iterations[i]
        ζ₃_max = maximum([ζ₃_max, maximum(ζ₃_xz[])])
    end

    ζ₃_max = minimum([ζ₃_max, 20*f])

    @info "Drawing first frame"

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    iter[] = 0
    fig = Figure(size = (1280, 720))
    ax_b = Axis(fig[1, 1][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Buoyancy, }b")
    ax_ζ = Axis(fig[1, 2][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Vertical vorticity, }\zeta/f")
    hm_b = heatmap!(ax_b, xb/1kilometer, zb/1kilometer, b_xz; colorrange = (-2.5*b_max, 1.5*b_max));
    hm_ζ₃ = heatmap!(ax_ζ, xζ₃/1kilometer, zζ₃/1kilometer, ζ_on_f; colormap = :coolwarm, colorrange = (-ζ₃_max/f, ζ₃_max/f));
    Colorbar(fig[1, 1][1, 2], hm_b)
    Colorbar(fig[1, 2][1, 2], hm_ζ₃)

    ##############################################
    # SEE IF CAN CHANGE COLORS, ADD TIME COUNTER #
    ##############################################

    display(fig)
    
    @info "Making an animation from saved data..."
    record(i -> iter[] = i,
           fig,
           "pretty_things/" * label * "_xz.mp4",
           iterations[Int64(round(length(iterations)*0.5)) : length(iterations)],
           framerate = 20)

end

nothing
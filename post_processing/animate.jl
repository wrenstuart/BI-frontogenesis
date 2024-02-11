using Oceananigans, JLD2, Plots, Printf
using Oceananigans.Units



function BI_plot(label)

    # Set the two dimensional parameters
    H = 50    # Depth of mixed layer
    f = 1e-4  # Coriolis parameter
    N² = 0

    filename_xy_top = "raw_data/" * label * "_BI_xy"

    # Read in the first iteration.  We do this to load the grid
    b_ic = FieldTimeSeries(filename_xy_top * ".jld2", "b", iterations = 0)
    ζ₃_ic = FieldTimeSeries(filename_xy_top * ".jld2", "ζ₃", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xb, yb, zb = nodes(b_ic)
    xζ₃, yζ₃, zζ₃ = nodes(ζ₃_ic)

    # Now, open the file with our data
    file_xy = jldopen(filename_xy_top * ".jld2")

    # Extract a vector of iterations
    iterations = parse.(Int, keys(file_xy["timeseries/t"]))
    # iterations[i] is the number of timesteps passed for the i-th frame

    @info "Making an animation from saved data..."

    t_save = zeros(length(iterations))  # Contains the actual time

    # Calculate the maximum relative vorticity and buoyancy flux to set the scale for the colourmap
    ζ₃_max = 0
    b_max = maximum(b_ic)
    b_min = 0

    for i = 11 : length(iterations)
        iter = iterations[i]
        t = file_xy["timeseries/t/$iter"]
        t_save[i] = t
        ζ₃_xy = file_xy["timeseries/ζ₃/$iter"][:, 1, :]
        this_ζ₃_max = maximum(ζ₃_xy)
        ζ₃_max = maximum([ζ₃_max, this_ζ₃_max])
    end

    # Here, we loop over all iterations
    anim = @animate for (i, iter) in enumerate(iterations[Int64(round(length(iterations)*0.4)) : length(iterations)])

        #@info "Drawing frame $i from iteration $iter..."

        b_xy = file_xy["timeseries/b/$iter"][:, :, 1];
        ζ₃_xy = file_xy["timeseries/ζ₃/$iter"][:, :, 1];

        t = file_xy["timeseries/t/$iter"];  # The time that has elapsed by this iteration

        # Save some variables to plot at the end
        t_save[i] = t # save the time

            b_xy_plot = Plots.heatmap(xb/1kilometer, yb/1kilometer, b_xy'/b_max; color = :viridis, xlabel = "\$x/\\mathrm{km}\$", ylabel = "\$y/\\mathrm{km}\$", clims=(-0.5, 1.5), ticklabel = []);
            ζ₃_xy_plot = Plots.heatmap(xζ₃/1kilometer, yζ₃/1kilometer, ζ₃_xy'/f; color = :balance, xlabel = "\$x/\\mathrm{km}\$", ylabel = "\$y/\\mathrm{km}\$", clims=(-ζ₃_max/f, ζ₃_max/f));

        b_title = @sprintf("\$b\$");
        ζ₃_title = @sprintf("\$ζ/f\$");

        # Combine the sub-plots into a single figure
        Plots.plot(b_xy_plot, ζ₃_xy_plot,
        layout = (1, 2),
        title = [b_title ζ₃_title], size = (1280, 720))#aspect_ratio = :equal)

        iter == iterations[end] && close(file_xy)
    end

    close(file_xy)

    # Save the animation to a file
    mp4(anim, "pretty_things/" * label * ".mp4", fps = 20) # hide

end
using CairoMakie
using Oceananigans
using JLD2
using Printf
using Oceananigans.Units
using Makie.Colors
using Statistics
using OffsetArrays
using FFTW

#=function snap_ζ_xy(label, frac)

    # Set the two dimensional parameters
    H = 50    # Depth of mixed layer
    f = 1e-4  # Coriolis parameter

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"

    # Read in the first iteration. We do this to load the grid
    ζ_ic = FieldTimeSeries(filename_xy_top, "ζ₃", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    x, y, z = nodes(ζ_ic)

    # Now, open the file with our data
    file = jldopen(filename_xy_top)

    iterations = parse.(Int, keys(file["timeseries/t"]))
    iterations_animated = iterations[Int64(round(length(iterations)*0.5)) : length(iterations)]
    iter = iterations_animated[Int(round(frac*length(iterations)))]

    # Set up observables for plotting that will update as the iteration number is updated
    ζ_xy = file["timeseries/ζ₃/$iter"][:,:,1]

    # Calculate the maximum relative vorticity and buoyancy flux to set the scale for the colourmap
    ζ_max = maximum(ζ_xy)
    ζ_max = minimum([ζ_max, 20f])

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    fig = Figure(resolution = (440, 440))
    #xtitle = Axis((fig[1,1]), ylabel = L"$y/\mathrm{km}$", yticks = ([], []))
    ax_ζ = Axis(fig[1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\zeta/f", aspect = 1)
    hm_ζ = heatmap!(ax_ζ, x/1kilometer, y/1kilometer, ζ_xy/f; colormap = :coolwarm, colorrange = (-ζ_max/f, ζ_max/f));
    Colorbar(fig[1, 2], hm_ζ)

    resize_to_layout!(fig)
    
    ##############################################
    # SEE IF CAN CHANGE COLORS, ADD TIME COUNTER #
    ##############################################

    display(fig)
    save("pretty_things/" * "ζ-snap.pdf", fig)

end

function snap_b_xy(label, frac)

    # Set the two dimensional parameters
    H = 50    # Depth of mixed layer
    f = 1e-4  # Coriolis parameter

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"

    # Read in the first iteration. We do this to load the grid
    ζ_ic = FieldTimeSeries(filename_xy_top, "b", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    x, y, z = nodes(ζ_ic)

    # Now, open the file with our data
    file = jldopen(filename_xy_top)

    iterations = parse.(Int, keys(file["timeseries/t"]))
    iterations_animated = iterations[Int64(round(length(iterations)*0.5)) : length(iterations)]
    iter = iterations_animated[Int(round(frac*length(iterations)))]

    # Set up observables for plotting that will update as the iteration number is updated
    ζ_xy = file["timeseries/b/$iter"][:,:,1]

    # Calculate the maximum relative vorticity and buoyancy flux to set the scale for the colourmap
    ζ_max = maximum(abs.(ζ_xy))

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    fig = Figure(resolution = (440, 440))
    #xtitle = Axis((fig[1,1]), ylabel = L"$y/\mathrm{km}$", yticks = ([], []))
    ax_ζ = Axis(fig[1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"b", aspect = 1)
    hm_ζ = heatmap!(ax_ζ, x/1kilometer, y/1kilometer, ζ_xy; colorrange = (-ζ_max, ζ_max));
    Colorbar(fig[1, 2], hm_ζ)

    resize_to_layout!(fig)
    
    ##############################################
    # SEE IF CAN CHANGE COLORS, ADD TIME COUNTER #
    ##############################################

    display(fig)
    save("pretty_things/" * "b-snap.pdf", fig)

end

function snap_δ_xy(label, frac)

    # Set the two dimensional parameters
    H = 50    # Depth of mixed layer
    f = 1e-4  # Coriolis parameter

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"

    # Read in the first iteration. We do this to load the grid
    b_ic = FieldTimeSeries(filename_xy_top, "δ", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    x, y, z = nodes(b_ic)

    # Now, open the file with our data
    file = jldopen(filename_xy_top)

    iterations = parse.(Int, keys(file["timeseries/t"]))
    iterations_animated = iterations[Int64(round(length(iterations)*0.5)) : length(iterations)]
    iter = iterations_animated[Int(round(frac*length(iterations)))]

    # Set up observables for plotting that will update as the iteration number is updated
    ζ_xy = file["timeseries/δ/$iter"][:,:,1]

    # Calculate the maximum relative vorticity and buoyancy flux to set the scale for the colourmap
    ζ_max = maximum(abs.(ζ_xy))
    ζ_max = minimum([ζ_max, 20f])

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    fig = Figure(resolution = (440, 440))
    #xtitle = Axis((fig[1,1]), ylabel = L"$y/\mathrm{km}$", yticks = ([], []))
    ax_ζ = Axis(fig[1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\delta/f", aspect = 1)
    hm_ζ = heatmap!(ax_ζ, x/1kilometer, y/1kilometer, colormap = :coolwarm, ζ_xy/f; colorrange = (-ζ_max/f, ζ_max/f));
    Colorbar(fig[1, 2], hm_ζ)

    resize_to_layout!(fig)
    
    ##############################################
    # SEE IF CAN CHANGE COLORS, ADD TIME COUNTER #
    ##############################################

    display(fig)
    save("pretty_things/" * "δ-snap.pdf", fig)

end=#

function snap_pic(label, ft, b_range, ζ_on_f_range, δ_on_f_range; x_labels = false, name = "none")

    # Set the two dimensional parameters
    H = 50    # Depth of mixed layer
    f = 1e-4  # Coriolis parameter

    filename = "raw_data/" * label * "_BI_xy" * ".jld2"

    # Read in the first iteration. We do this to load the grid
    ζ_ic = FieldTimeSeries(filename, "ζ", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    x, y, ~ = nodes(ζ_ic)

    # Now, open the file with our data
    file = jldopen(filename)

    iterations = parse.(Int, keys(file["timeseries/t"]))
    ts = [file["timeseries/t/$iter"] for iter in iterations]
    iter = iterations[argmin(abs.(ts .- ft/f))]

    # Set up observables for plotting that will update as the iteration number is updated
    b = file["timeseries/b/$iter"][:,:,1]
    ζ = file["timeseries/ζ/$iter"][:,:,1]
    δ = file["timeseries/δ/$iter"][:,:,1]
    filename_drifters = "raw_data/" * label * "_particles.jld2"
    file_drifters = jldopen(filename_drifters)
    drifters_x = file_drifters["timeseries/particles/$iter"].x
    drifters_y = file_drifters["timeseries/particles/$iter"].y

    fig = Figure()
    if x_labels
        ax_b = Axis(fig[1, 1][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Buoyancy, }b", width = 250, height = 250)
        ax_ζ = Axis(fig[1, 2][1, 1], xlabel = L"$x/\mathrm{km}$", title = L"\text{Vertical vorticity, }\zeta/f", width = 250, height = 250)
        ax_δ = Axis(fig[1, 3][1, 1], xlabel = L"$x/\mathrm{km}$", title = L"\text{Horizontal divergence, }\delta/f", width = 250, height = 250)
    else
        ax_b = Axis(fig[1, 1][1, 1], ylabel = L"$y/\mathrm{km}$", title = L"\text{Buoyancy, }b", width = 250, height = 250)
        ax_ζ = Axis(fig[1, 2][1, 1], title = L"\text{Vertical vorticity, }\zeta/f", width = 250, height = 250)
        ax_δ = Axis(fig[1, 3][1, 1], title = L"\text{Horizontal divergence, }\delta/f", width = 250, height = 250)
    end
    hm_b = heatmap!(ax_b, x/1kilometer, y/1kilometer, b; colorrange = b_range)
    hm_ζ = heatmap!(ax_ζ, x/1kilometer, y/1kilometer, ζ/f; colormap = :coolwarm, colorrange = ζ_on_f_range)
    hm_δ = heatmap!(ax_δ, x/1kilometer, y/1kilometer, δ/f; colormap = :coolwarm, colorrange = δ_on_f_range)
    scatter!(ax_b, drifters_x/1kilometer, drifters_y/1kilometer, marker = '.', markersize = 20, color = :black)
    scatter!(ax_ζ, drifters_x/1kilometer, drifters_y/1kilometer, marker = '.', markersize = 20, color = :black)
    scatter!(ax_δ, drifters_x/1kilometer, drifters_y/1kilometer, marker = '.', markersize = 20, color = :black)
    
    Colorbar(fig[1, 1][1, 2], hm_b)
    Colorbar(fig[1, 2][1, 2], hm_ζ)
    Colorbar(fig[1, 3][1, 2], hm_δ)
    Makie.Label(fig[0, 1:3], L"ft=%$(ft)")

    resize_to_layout!(fig)
    display(fig)
    if name ≠ "none"
        save("pretty_things/" * name * ".png", fig)
    end

end

function gaussian_filter_2d(z, m_cut, n_cut)
    (M, N) = size(z)
    z_f = fft(z)
    for i in 1:M, j in 1:N
        m = minimum([mod(i-1,M), mod(1-i,M)])
        n = minimum([mod(j-1,N), mod(1-j,N)])
        z_f[i, j] *= exp(-(m^2/m_cut^2 + n^2/n_cut^2))
    end
    real(ifft(z_f))
end

function per(i,N)
    # For periodic arrays
    mod(i-1, N) + 1
end

function snap_fdetect(label, ft, ∇b_scale = 5e-6, L_scale = 8000)

    f = 1e-4
    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"

    # Read in the first iteration. We do this to load the grid
    b_ic = FieldTimeSeries(filename_xy_top, "b", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xb, yb, ~ = nodes(b_ic)
    M = length(xb)
    N = length(yb)
    Δx = xb[2] - xb[1]
    Δy = yb[2] - yb[1]
    m_cut = Int(round((xb[2] - xb[1]) * M / L_scale))
    n_cut = Int(round((yb[2] - yb[1]) * N / L_scale))

    # Now, open the file with our data
    file = jldopen(filename_xy_top)
    iterations = parse.(Int, keys(file["timeseries/t"]))
    ts = [file["timeseries/t/$iter"] for iter in iterations]
    iter = iterations[argmin(abs.(ts .- ft/f))]

    frames = 1:length(iterations)
    front_highlight = OffsetArrays.no_offset_view(zeros(frames, M, N))
    front_diagnose = OffsetArrays.no_offset_view(zeros(frames, M, N))

    b = file["timeseries/b/$iter"][:, :, 1]
    b_x = file["timeseries/b_x/$iter"][:, :, 1]
    b_y = file["timeseries/b_y/$iter"][:, :, 1]
    abs∇b = [(x < 1e-4 ? x : 0) for x in (b_x.^2 + b_y.^2) .^ 0.5]
    #∇²b = [(b_x[per(i+1,M),j] - b_x[per(i-1,M),j])/2Δx + (b_y[i,per(j+1,N)] - b_y[i,per(j-1,N)])/2Δy for i in 1:M, j in 1:N]
    ∇b_filter = [(x > ∇b_scale ? 1 : 0) for x in abs∇b]
    #∇²b_filter = [(x > 10 * ∇b_scale/L_scale ? 1 : 0) for x in ∇²b]
    front_filter = ∇b_filter# .| ∇²b_filter
    front_highlight = front_filter

    filt_abs∇b = gaussian_filter_2d(abs∇b, m_cut, n_cut)
    b_x_filt = gaussian_filter_2d(b_x, m_cut, n_cut)
    b_y_filt = gaussian_filter_2d(b_y, m_cut, n_cut)
    abs_filt∇b = (b_x_filt.^2 + b_y_filt.^2) .^ 0.5
    𝒻 = abs_filt∇b ./ filt_abs∇b
    front_diagnose = 𝒻 .* front_filter

    fig = Figure()
    ax_b = Axis(fig[1, 1][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Buoyancy, }b", width = 250, height = 250)
    hm_b = heatmap!(ax_b, xb/1kilometer, yb/1kilometer, b; colorrange = (minimum(b), maximum(b)))
    ax_∇b = Axis(fig[1, 2][1, 1], xlabel = L"$x/\mathrm{km}$", title = L"\mathcal{F}", width = 250, height = 250)
    hm_∇b = heatmap!(ax_∇b, xb/1kilometer, yb/1kilometer, front_diagnose; colorrange = (0, 1), colormap = :thermal)
    Colorbar(fig[1, 1][1, 2], hm_b, height = 250)
    Colorbar(fig[1, 2][1, 2], hm_∇b, height = 250)

    resize_to_layout!(fig)
    display(fig)
    save("pretty_things/" * "fdetect-snap.png", fig)
    
end
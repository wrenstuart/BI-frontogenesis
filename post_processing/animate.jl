using CairoMakie
using Oceananigans
using JLD2
using Printf
using Oceananigans.Units
using Makie.Colors
using Statistics
using OffsetArrays
using FFTW

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

    δ = lift(iter -> file["timeseries/δ/$iter"][:, :, 1], iter)
    u_x = lift(iter -> file["timeseries/u_x/$iter"][:, :, 1], iter)
    v_x = lift(iter -> file["timeseries/v_x/$iter"][:, :, 1], iter)
    u_y = lift((v_x, ζ) -> v_x - ζ, v_x, ζ₃_xy)
    v_y = lift((u_x, δ) -> δ - u_x, u_x, δ)
    e = lift((u_x, u_y, v_x, v_y) -> u_x .* v_y - (u_y + v_x).^2 / 4, u_x, u_y, v_x, v_y)

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
    #hm_ζ₃ = heatmap!(ax_ζ, xζ₃/1kilometer, yζ₃/1kilometer, e; colormap = :coolwarm, colorrange = (-5f^2, 5f^2));
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

function ani_zeta_hist(label)

    # Set the two dimensional parameters
    H = 50    # Depth of mixed layer
    f = 1e-4  # Coriolis parameter

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"

    # Now, open the file with our data
    file = jldopen(filename_xy_top)

    # Set up observables for plotting that will update as the iteration number is updated
    iter = Observable(0)
    ζ₃_xy = lift(iter -> file["timeseries/ζ₃/$iter"][:, :, 1], iter) # Surface vertical vorticity at this iteration
    ζ_on_fs = lift(iter -> vec(ζ₃_xy[])/f, iter)
    t = lift(iter -> file["timeseries/t/$iter"], iter)   # Time elapsed by this iteration

    # Extract the values that iter can take
    iterations = parse.(Int, keys(file["timeseries/t"]))

    @info "Drawing first frame"

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    iter[] = 0
    fig = Figure(size = (1280, 720))
    ax = Axis(fig[1, 1][1, 1], xlabel = L"$\zeta/f$", ylabel = L"\text{Frequency density}", limits = ((-5, 10), (0, 2)))
    h = hist!(ax, ζ_on_fs, bins = -5:0.025:15, normalization = :pdf)

    ##############################################
    # SEE IF CAN CHANGE COLORS, ADD TIME COUNTER #
    ##############################################

    display(fig)
    
    @info "Making an animation from saved data..."
    record(i -> iter[] = i,
           fig,
           "pretty_things/" * label * "_hist.mp4",
           iterations[Int64(round(length(iterations)*0.5)) : length(iterations)],
           framerate = 20)

end

function ani_zeta_hist_cf(label₁, label₂)

    # Set the two dimensional parameters
    H = 50    # Depth of mixed layer
    f = 1e-4  # Coriolis parameter

    filename₁_xy_top = "raw_data/" * label₁ * "_BI_xy" * ".jld2"
    filename₂_xy_top = "raw_data/" * label₂ * "_BI_xy" * ".jld2"

    # Now, open the file with our data
    file₁ = jldopen(filename₁_xy_top)
    file₂ = jldopen(filename₂_xy_top)

    # Extract the values that iter can take
    iterations₁ = parse.(Int, keys(file₁["timeseries/t"]))
    iterations₂ = parse.(Int, keys(file₂["timeseries/t"]))
    t₁s = [file₁["timeseries/t/$iter"] for iter in iterations₁]
    t₂s = [file₂["timeseries/t/$iter"] for iter in iterations₂]

    ζ₁_scale = 0
    for iter in iterations₁[Int64(round(length(iterations₁)*0.5)) : length(iterations₁)]
        ζ₁_scale = maximum([mean(file₁["timeseries/ζ₃/$iter"][:, :, 1] .^ 2) ^ 0.5, ζ₁_scale])
    end

    t₁_ref = 0
    for iter in iterations₁[Int64(round(length(iterations₁)*0.5)) : length(iterations₁)]
        if mean(file₁["timeseries/ζ₃/$iter"][:, :, 1] .^ 2) >= ζ₁_scale^2 / 2
            t₁_ref = file₁["timeseries/t/$iter"]
            break
        end
    end
    t₂_ref = 0
    for iter in iterations₂[Int64(round(length(iterations₂)*0.5)) : length(iterations₂)]
        if mean(file₂["timeseries/ζ₃/$iter"][:, :, 1] .^ 2) >= ζ₁_scale^2 / 2
            t₂_ref = file₂["timeseries/t/$iter"]
            break
        end
    end
    Δt = t₂_ref - t₁_ref
    @info Δt
    @info Δt/t₁s[end]

    # Set up observables for plotting that will update as the iteration number is updated
    iter = Observable(0)
    t = lift(iter -> file₁["timeseries/t/$iter"], iter)   # Time elapsed by this iteration
    ζ₁ = lift(iter -> file₁["timeseries/ζ₃/$iter"][:, :, 1], iter) # Surface vertical vorticity at this iteration
    ζ₁_on_fs = lift(ζ₁ -> vec(ζ₁)/f, ζ₁)
    iter₂ = lift(iter -> iterations₂[findmin(abs.(t₂s .- (t[] + Δt)))[2]], iter)
    ζ₂ = lift(iter -> file₂["timeseries/ζ₃/$iter"][:, :, 1], iter₂) # Surface vertical vorticity at this iteration
    ζ₂_on_fs = lift(ζ₂ -> vec(ζ₂)/f, ζ₂)

    @info "Drawing first frame"

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    iter[] = 0
    fig = Figure(size = (1280, 720))
    ax = Axis(fig[1, 1][1, 1], xlabel = L"$\zeta/f$", ylabel = L"\text{Frequency density}", limits = ((-5, 10), (0, 2)))
    h = hist!(ax, ζ₁_on_fs, bins = -5:0.025:15, normalization = :pdf)
    h = hist!(ax, ζ₂_on_fs, bins = -5:0.025:15, normalization = :pdf)

    ##############################################
    # SEE IF CAN CHANGE COLORS, ADD TIME COUNTER #
    ##############################################

    display(fig)

    @info "Making an animation from saved data..."
    record(i -> iter[] = i,
           fig,
           "pretty_things/" * label₁ * label₂ * "_histcf.mp4",
           iterations₁[Int64(round(length(iterations₁)*0.5)) : length(iterations₁)],
           framerate = 20)

end

function low_pass_filter_2d(z, m_cut, n_cut)
    (M, N) = size(z)
    z_f = fft(z)
    z_f[1+m_cut:M-m_cut, :] .= 0
    z_f[:, 1+n_cut:N-n_cut] .= 0
    real(ifft(z_f))
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

function front_detection(label, ∇b_scale = 5e-6, L_scale = 8000)

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"

    # Read in the first iteration. We do this to load the grid
    b_ic = FieldTimeSeries(filename_xy_top, "b", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xb, yb, zb = nodes(b_ic)
    M = length(xb)
    N = length(yb)
    Δx = xb[2] - xb[1]
    Δy = yb[2] - yb[1]
    m_cut = Int(round((xb[2] - xb[1]) * M / L_scale))
    n_cut = Int(round((yb[2] - yb[1]) * N / L_scale))

    # Now, open the file with our data
    file = jldopen(filename_xy_top)
    iterations_full = parse.(Int, keys(file["timeseries/t"]))
    iterations = parse.(Int, keys(file["timeseries/t"]))[Int64(round(length(iterations_full)*0.5)) : length(iterations_full)]
    frames = 1:length(iterations)
    front_highlight = OffsetArrays.no_offset_view(zeros(frames, M, N))
    front_diagnose = OffsetArrays.no_offset_view(zeros(frames, M, N))

    for frame in frames

        @info string(frame) * "/" * string(length(frames))

        iter = iterations[frame]
        b_x = file["timeseries/b_x/$iter"][:, :, 1]
        b_y = file["timeseries/b_y/$iter"][:, :, 1]
        abs∇b = [(x < 1e-4 ? x : 0) for x in (b_x.^2 + b_y.^2) .^ 0.5]
        #∇²b = [(b_x[per(i+1,M),j] - b_x[per(i-1,M),j])/2Δx + (b_y[i,per(j+1,N)] - b_y[i,per(j-1,N)])/2Δy for i in 1:M, j in 1:N]
        ∇b_filter = [(x > ∇b_scale ? 0 : 0) for x in abs∇b]
        #∇²b_filter = [(x > 10 * ∇b_scale/L_scale ? 1 : 0) for x in ∇²b]
        front_filter = ∇b_filter# .| ∇²b_filter
        front_highlight[frame, :, :] = front_filter

        filt_abs∇b = gaussian_filter_2d(abs∇b, m_cut, n_cut)
        b_x_filt = gaussian_filter_2d(b_x, m_cut, n_cut)
        b_y_filt = gaussian_filter_2d(b_y, m_cut, n_cut)
        abs_filt∇b = (b_x_filt.^2 + b_y_filt.^2) .^ 0.5
        𝒻 = abs_filt∇b ./ filt_abs∇b
        front_diagnose[frame, :, :] = 𝒻 .* front_filter

    end

    frame = Observable(1)
    this_front_diagnose = lift(frame -> front_diagnose[frame, :, :], frame)
    this_front_highlight = lift(frame -> front_highlight[frame, :, :], frame)
    fig = Figure()
    ax_∇b = Axis(fig[1, 1][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\nabla b\text{ detection}")
    hm_b = heatmap!(ax_∇b, xb/1kilometer, yb/1kilometer, this_front_diagnose; colorrange = (0, 1));

    display(fig)
    
    @info "Making an animation from saved data..."
    record(i -> frame[] = i,
           fig,
           "pretty_things/" * label * "_fdetect" * ".mp4",
           frames,
           framerate = 20)
    
end

#####################################################################
# TRY REWRITING WITH OBSERVABLES ONLY MENTIONED IN THE PLOTTING BIT #
#####################################################################

nothing
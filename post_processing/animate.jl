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

function low_pass_filter_2d(z, m, n)
    z_f = fft(z)
    z_f[1+m_cut:M-m_cut, :] .= 0
    z_f[:, 1+n_cut:N-n_cut] .= 0
    z_smooth = real(ifft(z_f))
end

function front_detection(label)

    filename_xy_top = "raw_data/" * label * "_BI_xy" * ".jld2"

    # Read in the first iteration. We do this to load the grid
    b_ic = FieldTimeSeries(filename_xy_top, "b", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xb, yb, zb = nodes(b_ic)

    # Now, open the file with our data
    file = jldopen(filename_xy_top)

    # Set up observables for plotting that will update as the iteration number is updated
    iter = Observable(0)
    b = lift(iter -> file["timeseries/b/$iter"][:, :, 1], iter)   # Surface buoyancy at this iteration
    b_x = lift(iter -> file["timeseries/b_x/$iter"][:, :, 1], iter)   # Surface buoyancy at this iteration
    b_y = lift(iter -> file["timeseries/b_y/$iter"][:, :, 1], iter)   # Surface buoyancy at this iteration
    𝒻 = lift(iter -> file["timeseries/b/$iter"][:, :, 1], iter)   # Surface buoyancy at this iteration
    t = lift(iter -> file["timeseries/t/$iter"], iter)   # Time elapsed by this iteration
    # Extract the values that iter can take
    iterations = parse.(Int, keys(file["timeseries/t"]))
    b_max = maximum(b_ic)

    @info "Drawing first frame"

    abs_∇b = lift((b_x, b_y) -> (b_x.^2 + b_y.^2).^0.5, b_x, b_y)
    vec_abs_∇b = lift(x -> vec(x), abs_∇b)
    alt_abs_∇b = lift((b_x, b_y) -> (b_x.^2 + (b_y.-mean(b_y)).^2).^0.5, b_x, b_y)
    clean_abs_∇b = lift(abs_∇b -> [(x < 1e-4 ? x : mean(abs_∇b)) for x in abs_∇b], abs_∇b)
    vec_clean_abs_∇b = lift(x -> vec(x), clean_abs_∇b)
    #############################################################################################
    # THE CLEANING NEEDS TO BE DONE DUE TO INSTABILITIES; HOPEFULLY NOT PRESENT IN FINAL THING! #
    #############################################################################################
    #=∇b_scale = lift(clean_abs_∇b -> begin
                                  μ = mean(clean_abs_∇b)
                                  σ = var(clean_abs_∇b) .^ 0.5
                                  #@info(sum([(x > μ+3σ ? 1 : 0) for x in clean_abs_∇b])/length(clean_abs_∇b))
                                  σ
                              end,
                              clean_abs_∇b)=#
    ∇b_scale = lift(x -> var(x)^0.5, vec_clean_abs_∇b)
    ΔB = b[][1,end] - b[][1,1]
    vec_Δb₁ = lift(b -> vec([b[i,j] - j/size(b)[2] * ΔB for i in 1:size(b)[1], j in 1:size(b)[2]]), b)
    vec_Δb₂ = lift(b -> vec([b[i,j] .- mean(b[:,j]) for i in 1:size(b)[1], j in 1:size(b)[2]]), b)
    ∇²b = lift(iter -> OffsetArrays.no_offset_view(Field(∂x(FieldTimeSeries(filename_xy_top, "b_x", iterations = 0)) + ∂y(FieldTimeSeries(filename_xy_top, "b_y", iterations = 0)))[:,:,end]), iter)
    Δb_scale = lift(x -> var(x)^0.5, vec_Δb₂)
    
    #=L_scale = lift((Δb, ∇b) -> Δb/∇b, Δb_scale, ∇b_scale)
    on(L_scale) do L
        print(L)
        print('\n')
    end=#
    L_scale = 2000
    m_cut = Int(round((xb[end] - xb[1])/L_scale))
    n_cut = Int(round((yb[end] - yb[1])/L_scale))
    filter(z) = low_pass_filter_2d(z, m_cut, n_cut)
    abs_filt_∇b = lift(b_x, b_y) do b_x, b_y
        b_y = [(x < 1e-4 ? x : mean(b_y)) for x in b_y]
        b_x_filt = filter(b_x)
        b_y_filt = filter(b_y)
        (b_x.^2 + b_y.^2) .^ 0.5
    end
    filt_abs_∇b = lift(∇b -> filter(∇b), clean_abs_∇b, clean_abs_∇b)
    front_highlight = lift(clean_abs_∇b) do mag_∇b
        μ = mean(mag_∇b)
        σ = var(mag_∇b) ^ 0.5
        [(abs(x - μ) > 3.5σ ? 1 : 0) for x in mag_∇b]
    end
    fig = Figure()
    #ax = Axis(fig[1, 1][1, 1])
    #h = hist!(ax, vec_abs_∇b, bins = 0:1e-7:1e-5, normalisation = :pdf)
    #h = hist!(ax, vec_Δb₂, bins = -0.05:0.001:0.05, normalisation = :pdf)
    ax_∇b = Axis(fig[1, 1][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\nabla b\text{ detection}")
    hm_b = heatmap!(ax_∇b, xb/1kilometer, yb/1kilometer, front_highlight; colorrange = (0, 1));

    display(fig)
    
    @info "Making an animation from saved data..."
    record(i -> iter[] = i,
           fig,
           "pretty_things/" * label * "_fdetect" * ".mp4",
           iterations[Int64(round(length(iterations)*0.5)) : length(iterations)],
           framerate = 20)
    
end

nothing
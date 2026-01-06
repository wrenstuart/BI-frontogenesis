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

function ft_display(t::AbstractFloat)
    ft = @sprintf("%.2f", f*t)
    return L"ft=%$(ft)"
end

function ani_xy(label::String, a::AbstractFloat, b::AbstractFloat)  # Animate b, ζ and δ at surface

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
    iter = Observable(0)
    iter_index = lift(iter -> findfirst(x -> x == iter, iterations), iter)
    ζ_on_f = lift(iter -> file["timeseries/ζ/$iter"][:, :, 1]/f, iter)
    δ_on_f = lift(iter -> file["timeseries/δ/$iter"][:, :, 1]/f, iter)
    b_xy = lift(iter -> file["timeseries/b/$iter"][:, :, 1], iter)
    t = lift(iter -> file["timeseries/t/$iter"], iter)   # Time elapsed by this iteration
    str_ft = lift(t -> ft_display(t), t)

    # Calculate the maximum relative vorticity and buoyancy to set the scale for the colourmap
    ζ_max = 0
    b_max = maximum(b_ic)
    for i = Int(round(length(iterations)/10)) : length(iterations)
        iter[] = iterations[i]
        ζ_max = maximum([ζ_max, maximum(ζ_on_f[]*f)])
    end
    ζ_max = minimum([ζ_max, 20f])

    iter[] = 0
    ~, drifters = extract_tracked_drifter_data(label)
    @info length(iterations)
    xs = lift(iter_index -> [drifters[i][iter_index].x for i in eachindex(drifters)], iter_index)
    ys = lift(iter_index -> [drifters[i][iter_index].y for i in eachindex(drifters)], iter_index)
    xs_dspl = lift(xs -> xs/1kilometer, xs)
    ys_dspl = lift(ys -> ys/1kilometer, ys)

    @info "Drawing first frame"

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    fig = Figure(size = (950, 320))
    ax_b = Axis(fig[1, 1][1, 1], xlabel = L"$x/\mathrm{km}$", ylabel = L"$y/\mathrm{km}$", title = L"\text{Buoyancy, }b", width = 200, height = 200)
    ax_ζ = Axis(fig[1, 2][1, 1], xlabel = L"$x/\mathrm{km}$", title = L"\text{Vertical vorticity, }\zeta/f", width = 200, height = 200)
    ax_δ = Axis(fig[1, 3][1, 1], xlabel = L"$x/\mathrm{km}$", title = L"\text{Horizontal divergence, }\delta/f", width = 200, height = 200)
    hm_b = heatmap!(ax_b, xb/1kilometer, yb/1kilometer, b_xy; colorrange = (-0.5*b_max, 1.5*b_max));
    sc_b = scatter!(ax_b, xs_dspl, ys_dspl, marker = '.', markersize = 15, color = :black, transparency = true)
    hm_ζ = heatmap!(ax_ζ, xζ/1kilometer, yζ/1kilometer, ζ_on_f; colormap = :coolwarm, colorrange = (-ζ_max/f, ζ_max/f))
    sc_ζ = scatter!(ax_ζ, xs_dspl, ys_dspl, marker = '.', markersize = 15, color = :black, transparency = true)
    hm_δ = heatmap!(ax_δ, xδ/1kilometer, yδ/1kilometer, δ_on_f; colormap = :coolwarm, colorrange = (-ζ_max/f, ζ_max/f))
    sc_δ = scatter!(ax_δ, xs_dspl, ys_dspl, marker = '.', markersize = 15, color = :black, transparency = true)
    Colorbar(fig[1, 1][1, 2], hm_b, height = 200)
    Colorbar(fig[1, 2][1, 2], hm_ζ, height = 200)
    Colorbar(fig[1, 3][1, 2], hm_δ, height = 200)
    Makie.Label(fig[0, 1:3], str_ft)

    display(fig)
    
    @info "Making an animation from saved data..."
    CairoMakie.record(i -> iter[] = i,
           fig,
           pp_dir(label) * "bζδ-top-vid.mp4",
           iterations[maximum([Int64(round(length(iterations) * a )), 1]) : 2 : Int64(round(length(iterations) * b))],
           framerate = 20)

end

ani_xy(label::String) = ani_xy(label::String, 0.0, 1.0)

function ani_zeta_hist(label::String)

    check_pp_dir(label)
    filename_xy_top = data_dir(label) * "BI_xy.jld2"

    # Now, open the file with our data
    file = jldopen(filename_xy_top)

    # Set up observables for plotting that will update as the iteration number is updated
    iter = Observable(0)
    ζ_xy = lift(iter -> file["timeseries/ζ/$iter"][:, :, 1], iter) # Surface vertical vorticity at this iteration
    ζ_on_fs = lift(iter -> vec(ζ_xy[])/f, iter)
    t = lift(iter -> file["timeseries/t/$iter"], iter)   # Time elapsed by this iteration
    str_ft = lift(t -> ft_display(t), t)

    # Extract the values that iter can take
    iterations = parse.(Int, keys(file["timeseries/t"]))

    @info "Drawing first frame"

    # Create the plot, starting at t = 0
    # This will be updated as the observable iter is updated
    iter[] = 0
    fig = Figure(size = (400, 400))
    ax = Axis(fig[1, 1], xlabel = L"$\zeta/f$", ylabel = L"\text{Frequency density}", limits = ((-5, 10), (0, 2)), width = 320)
    hist!(ax, ζ_on_fs, bins = -5:0.025:15, normalization = :pdf)
    Makie.Label(fig[0, 1], str_ft)

    display(fig)
    
    @info "Making an animation from saved data..."
    record(i -> iter[] = i,
           fig,
           pp_dir(label) * "ζ-top-hist-vid.mp4",
           iterations[Int64(round(length(iterations)*0.5)) : length(iterations)],
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

function front_detection(label, ∇b_scale = 5e-6, L_scale = 8000)

    check_pp_dir(label)
    filename_xy_top = data_dir(label) * "BI_xy.jld2"

    # Read in the first iteration. We do this to load the grid
    b_ic = FieldTimeSeries(filename_xy_top, "b", iterations = 0)

    # Load in co-ordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xb, yb, ~ = nodes(b_ic)

    M = length(xb)
    N = length(yb)
    Δx = xb[2] - xb[1]
    Δy = yb[2] - yb[1]
    m_cut = Int(round(Δx * M / L_scale))    # Cutoff x-wavenumber (ignoring factors of π)
    n_cut = Int(round(Δy * N / L_scale))    # Cutoff y-wavenumber (ignoring factors of π)

    # Now, open the file with our data
    file = jldopen(filename_xy_top)
    iterations_full = parse.(Int, keys(file["timeseries/t"]))
    iterations = parse.(Int, keys(file["timeseries/t"]))[Int64(round(length(iterations_full)*0.5)) : length(iterations_full)]
    frames = 1:length(iterations)
    front_highlight = OffsetArrays.no_offset_view(zeros(frames, M, N))
    # Stores data as to whether something has been detected as a front of not
    front_diagnose = OffsetArrays.no_offset_view(zeros(frames, M, N))
    # Stores quantity that diagnoses front vs. filament

    for frame in frames

        @info string(frame) * "/" * string(length(frames))

        iter = iterations[frame]
        b_x = file["timeseries/b_x/$iter"][:, :, 1]
        b_y = file["timeseries/b_y/$iter"][:, :, 1]
        abs∇b = [(x < 1e-4 ? x : 0) for x in (b_x.^2 + b_y.^2) .^ 0.5]
        #∇²b = [(b_x[per(i+1,M),j] - b_x[per(i-1,M),j])/2Δx + (b_y[i,per(j+1,N)] - b_y[i,per(j-1,N)])/2Δy for i in 1:M, j in 1:N]
        ∇b_filter = [(x > ∇b_scale ? 0 : 1) for x in abs∇b]
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
           pp_dir(label) * "top-fdetect.mp4",
           frames,
           framerate = 20)
    
end

function ani_drifters(label::String)
    t_drifters, drifters = get_drifter_data(label)
    ani_drifters(label, drifters[1])
end

function ani_drifters(label::String, drifter)     # Quite different to pre-re-factored function of same name
    
    fig = Figure(size = (950, 950))

    data = topdata(label)
    iterations = parse.(Int, keys(data.file["timeseries/t"]))
    t_drifters, ~ = extract_tracked_drifter_data(label)
    ft = f * t_drifters

    frame = Observable(1)
    iter = lift(i -> iterations[i], frame)
    t_obs = lift(iter -> data.file["timeseries/t/$iter"], iter)
    ζ_on_f = lift(i -> data.file["timeseries/ζ/$i"][:, :, 1]/f, iter)
    δ_on_f = lift(i -> data.file["timeseries/δ/$i"][:, :, 1]/f, iter)
    b = lift(i -> data.file["timeseries/b/$i"][:, :, 1], iter)
    i_drifter = lift(t -> argmin(abs.(t_drifters .- t)), t_obs)
    ft_obs = lift(i -> ft[i], i_drifter)
    tracers_now_x = lift(i -> drifter[i].x/1e3, i_drifter)
    tracers_now_y = lift(i -> drifter[i].y/1e3, i_drifter)

    b_ic = data.file["timeseries/b/0"][:, :, 1]
    b_max = maximum(b_ic)

    ax_ζ = Axis(fig[1, 1][1, 1], aspect = 1)
    hm_ζ = heatmap!(ax_ζ, data.xᶠ/1e3, data.yᶠ/1e3, ζ_on_f, colormap = :coolwarm, colorrange = (-20, 20), height = 200);
    scatter!(ax_ζ, tracers_now_x, tracers_now_y, marker = '.', markersize = 15, color = :black)
    Colorbar(fig[1, 1][1, 2], hm_ζ, height = 200)
    ax_δ = Axis(fig[1, 2][1, 1], aspect = 1)
    hm_δ = heatmap!(ax_δ, data.xᶜ/1e3, data.yᶜ/1e3, δ_on_f, colormap = :coolwarm, colorrange = (-20, 20), height = 200);
    scatter!(ax_δ, tracers_now_x, tracers_now_y, marker = '.', markersize = 15, color = :black)
    Colorbar(fig[1, 2][1, 2], hm_δ, height = 200)
    ax_b = Axis(fig[1, 3][1, 1], aspect = 1)
    hm_b = heatmap!(ax_b, data.xᶜ/1e3, data.yᶜ/1e3, b, colorrange = (-0.5b_max, 1.5b_max), height = 200);
    scatter!(ax_b, tracers_now_x, tracers_now_y, marker = '.', markersize = 15, color = :black)
    Colorbar(fig[1, 3][1, 2], hm_b, height = 200)
    
    resize_to_layout!(fig)
    display(fig)

    CairoMakie.record(i -> frame[] = i, fig, "pretty_things/tracer_" * label * ".mp4", 1 : length(iterations), framerate = 20)
    
end

################################
# HERE ON YET TO BE REFACTORED #
################################

function ζ_δ_lagr(label)

    # Makes a scatter plot of ζ and δ tracked by drifters, animated with time

    data = topdata(label)
    drifters = extract_tracers(label)
    drifters = [drifter for drifter in drifters[1:4:end]]
    ζs_t = zeros((length(drifters), length(drifters[1])))
    @info size(ζs_t)
    δs_t = zeros((length(drifters), length(drifters[1])))
    # Each of the above is indexed by drifter number, then iteration number

    for (i, drifter) in enumerate(drifters)
        @info i
        ~, ζ_t = @time lagr_track(data, plotting_vars.ζ_on_f, drifter)
        ~, δ_t = @time lagr_track(data, plotting_vars.δ_on_f, drifter)
        ζs_t[i, :] = ζ_t
        δs_t[i, :] = δ_t
    end

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = 1, limits = ((-20, 20), (-20, 20)))
    frame = Observable(1)
    ζs = lift(i -> ζs_t[:, i], frame)
    δs = lift(i -> δs_t[:, i], frame)
    scatter!(ax, ζs, δs, marker = '.', markersize = 30, color = :black)

    record(i -> frame[] = i, fig, "pretty_things/ζ-δ-drifter_" * label * ".mp4", 1 : length(drifters[1]), framerate = 20)

end
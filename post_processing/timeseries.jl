include("pp-io.jl")
include("drifters-refactored.jl")

label = "test_extra_visc_low_res"

function investigate_lagr_ζ_balance(label::String, drifter_num::Int64, plot_mode = "interpolated")
    
    # WE ASSUME THAT BI_XY AND PARTICLE ITERATIONS ARE THE SAME

    check_pp_dir(label)
    eul_data = topdata(label)
    t, tracked_drifter_data = extract_tracked_drifter_data(label)
    num_iters = length(tracked_drifter_data[drifter_num])
    iterations = eul_data.iterations

    grid_pos = (Face(), Face())

    x = [tracked_drifter_data[drifter_num][i].x for i = 1 : num_iters]
    y = [tracked_drifter_data[drifter_num][i].y for i = 1 : num_iters]
    ζ_tracked = [tracked_drifter_data[drifter_num][i].ζ for i = 1 : num_iters]
    ζ_tendency_tracked = [tracked_drifter_data[drifter_num][i].ζ_tendency for i = 1 : num_iters]
    ζ_cor_tracked = [tracked_drifter_data[drifter_num][i].ζ_cor for i = 1 : num_iters]
    ζ_visc_tracked = [tracked_drifter_data[drifter_num][i].ζ_visc for i = 1 : num_iters]
    ζ_err_tracked = [tracked_drifter_data[drifter_num][i].ζ_err for i = 1 : num_iters]
    F_ζ_hor_tracked = [tracked_drifter_data[drifter_num][i].F_ζ_hor for i = 1 : num_iters]
    F_ζ_vrt_tracked = [tracked_drifter_data[drifter_num][i].F_ζ_vrt for i = 1 : num_iters]
    ζ_adv_tracked = [tracked_drifter_data[drifter_num][i].ζ_adv for i = 1 : num_iters]    
    ζ_interpolated = extract_interpolated_drifter_data(eul_data, "ζ", grid_pos, x, y, t)
    ζ_tendency_interpolated = extract_interpolated_drifter_data(eul_data, "ζ_tendency", grid_pos, x, y, t)
    ζ_cor_interpolated = extract_interpolated_drifter_data(eul_data, "ζ_cor", grid_pos, x, y, t)
    ζ_visc_interpolated = extract_interpolated_drifter_data(eul_data, "ζ_visc", grid_pos, x, y, t)
    ζ_err_interpolated = extract_interpolated_drifter_data(eul_data, "ζ_err", grid_pos, x, y, t)
    F_ζ_hor_interpolated = extract_interpolated_drifter_data(eul_data, "F_ζ_hor", grid_pos, x, y, t)
    F_ζ_vrt_interpolated = extract_interpolated_drifter_data(eul_data, "F_ζ_vrt", grid_pos, x, y, t)
    ζ_adv_interpolated = extract_interpolated_drifter_data(eul_data, "ζ_adv", grid_pos, x, y, t)

    fig = Figure()
    ax = Axis(fig[1, 1])
    if plot_mode == "tracked"
        lines!(f*t, ζ_tendency_tracked + ζ_adv_tracked + ζ_err_tracked, label = L"\mathrm{D}\zeta/\mathrm{D}t")
        lines!(f*t, ζ_cor_tracked, label = L"\zeta_\text{Cor}")
        lines!(f*t, ζ_visc_tracked, label = L"\zeta_\text{visc}")
        # lines!(f*t, ζ_err_tracked, label = L"\zeta_\text{err}")
        lines!(f*t, F_ζ_hor_tracked, label = L"F_{\zeta,\text{hor}}")
        # lines!(f*t, F_ζ_vrt_tracked, label = L"F_{\zeta,\text{vrt}}")
        lines!(f*t, ζ_tendency_tracked + ζ_adv_tracked - (
            ζ_cor_tracked + ζ_visc_tracked + ζ_err_tracked + F_ζ_hor_tracked),
            label = "residual")
    elseif plot_mode == "interpolated"
        lines!(f*t, ζ_tendency_interpolated + ζ_adv_interpolated + ζ_err_interpolated, label = L"\mathrm{D}\zeta/\mathrm{D}t")
        lines!(f*t, ζ_cor_interpolated, label = L"\zeta_\text{Cor}")
        lines!(f*t, ζ_visc_interpolated, label = L"\zeta_\text{visc}")
        # lines!(f*t, ζ_err_interpolated, label = L"\zeta_\text{err}")
        lines!(f*t, F_ζ_hor_interpolated, label = L"F_{\zeta,\text{hor}}")
        # lines!(f*t, F_ζ_vrt_interpolated, label = L"F_{\zeta,\text{vrt}}")
        lines!(f*t, ζ_tendency_interpolated + ζ_adv_interpolated + ζ_err_interpolated - (
            ζ_cor_interpolated + ζ_visc_interpolated + F_ζ_hor_interpolated),
            label = "residual")
    end
    ylims!(ax, -4e-7, 4e-7)
    axislegend()
    display(fig)

end

lagr_ζ_balance(label::String)  = lagr_ζ_balance(label, 1)
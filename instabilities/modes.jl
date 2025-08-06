using Printf
using Statistics
include("diff-ops.jl")

function amplitude(u, v, w, b, f, aspect)
    
    return (mean(abs.(u).^2) + mean(abs.(v).^2) + mean(abs.(w/aspect).^2) + mean(abs.(aspect/f * b).^2)) .^ 0.5

end

function least_stable_mode(Ri, k, l; H = 50, f = 1e-4, νᵥ = 1e-3, νₕ = 1e+1, N² = 1e-4, N = 50, rate_only = false)
    
    M² = (N²*f^2/Ri) ^ 0.5
    Δz = H/(N-1)
    z = vec(-H : Δz : 0)

    Id = Id_mat(N)
    O = 0Id
    D = diff_mat(1, N, Δz, error_ord = 2)
    D² = diff_mat(2, N, Δz, error_ord = 2)
    i = im
    κ² = k^2 + l^2
    U_vec = -M²/f * (z.+H)
    U = U_vec .* Id
    ikU = i*k*U
    ik = i*k*Id
    il = i*l*Id
    𝒟 = νᵥ*D² - νₕ*κ²*Id
    ℒ = 𝒟 - ikU

    𝒜 = [Id   O    O    O    O;
        O    Id   O    O    O;
        O    O    Id   O    O;
        O    O    O    Id   O;
        O    O    O    O    O]

    ℬ = [ℒ         f*Id      M²/f*Id   O        -ik;
        -f*Id      ℒ         O         O        -il;
        O         O         ℒ         Id       -D;
        O        -M²*Id    -N²*Id     ℒ         O;
        ik        il        D         O         O]

    Id₀ = Id[1:1, :]
    Id₁ = Id[N:N, :]
    D₀ = D[1:1, :]
    D₁ = D[N:N, :]
    O_ = O[1:1, :]
    𝒞 = [D₀  O_  O_  O_  O_;
        D₁  O_  O_  O_  O_;
        O_  D₀  O_  O_  O_;
        O_  D₁  O_  O_  O_;
        O_  O_  Id₀ O_  O_;
        O_  O_  Id₁ O_  O_;
        O_  O_  O_  D₀  O_;
        O_  O_  O_  D₁  O_]
    removed = [1, N, N+1, 2N, 2N+1, 3N, 3N+1, 4N]

    σ, 𝐯 = evals_with_constraints(𝒜, ℬ, 𝒞, removed)

    function number_bad(σ)
        n_cutoff = 0
        for σ in σ
            if isnan(σ)
                n_cutoff += 1
            elseif real(σ) > 1
                n_cutoff += 1
            end
        end
        return n_cutoff
    end
    n_cutoff = number_bad(σ)
    σ = σ[1:end-n_cutoff]
    𝐯 = 𝐯[:, 1:length(σ)]

    index = length(σ)

    if rate_only
        return σ[index]
    end

    u = 𝐯[0N+1:1N, index]
    v = 𝐯[1N+1:2N, index]
    w = 𝐯[2N+1:3N, index]
    b = 𝐯[3N+1:4N, index]
    p = 𝐯[4N+1:5N, index]

    A = amplitude(u, v, w, b, f, ((k^2+l^2)^0.5 * H)/2π)

    return u/A, v/A, w/A, b/A, σ[index]

end

function generate_ic(Ri, L, U; N = 50, H = 50)
    
    u_modes = []
    v_modes = []
    w_modes = []
    b_modes = []
    # ~, ~, ~, ~, σ₀ = least_stable_mode(Ri, 4π/L, 0, N = N)
    for i = 0:10, j = 0:10
        if !(i == 0 && j == 0)
            k = i * 2π/L
            l = j * 2π/L
            û, v̂, ŵ, b̂, σ = least_stable_mode(Ri, k, l, N = N)
            A = U * ((i==2 && j==0) ? 0.03 : 0.0005 * rand() * exp(2π*im*rand()))
            # A = U * 0.01 * 10000^(real(σ)/real(σ₀) - 1) * exp(2π*im*rand())
            # @info i, j, real(σ)/real(σ₀) - 1
            push!(u_modes, ((k, l), A * û))
            push!(v_modes, ((k, l), A * v̂))
            push!(w_modes, ((k, l), A * ŵ))
            push!(b_modes, ((k, l), A * b̂))
        end
    end

    function interpolate_mode(f, x)
        # Assuming x ∈ [0, 1] with f[1] = f(0) and f[end] = f(1)
        i = 1 + (length(f)-1) * x
        if isinteger(i)
            return f[Int(i)]
        else
            i₀ = Int(floor(i))
            i₁ = Int(ceil(i))
            return (i₁-i) * f[i₀] + (i-i₀) * f[i₁]
        end
    end

    function uᵢ(x, y, z)
        val = 0
        for mode in u_modes
            k = mode[1][1]
            l = mode[1][2]
            û = mode[2]
            val += real(exp(im*(k*x+l*y)) * interpolate_mode(û, z/H+1))
        end
        return val + 0.001 * U * randn()
    end
    function vᵢ(x, y, z)
        val = 0
        for mode in v_modes
            k = mode[1][1]
            l = mode[1][2]
            v̂ = mode[2]
            val += real(exp(im*(k*x+l*y)) * interpolate_mode(v̂, z/H+1))
        end
        return val + 0.001 * U * randn()
    end
    function wᵢ(x, y, z)
        val = 0
        for mode in w_modes
            k = mode[1][1]
            l = mode[1][2]
            ŵ = mode[2]
            val += real(exp(im*(k*x+l*y)) * interpolate_mode(ŵ, z/H+1))
        end
        return val + 0.001 * U * randn()
    end
    function bᵢ(x, y, z)
        val = 0
        for mode in b_modes
            k = mode[1][1]
            l = mode[1][2]
            b̂ = mode[2]
            val += real(exp(im*(k*x+l*y)) * interpolate_mode(b̂, z/H+1))
        end
        return val
    end

    return uᵢ, vᵢ, wᵢ, bᵢ

end

#=z = -50:0
u, v, w, b = generate_ic(1, 28000, 0.5, N = 50)
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, u.(0, 0, z), z)
lines!(ax, v.(0, 0, z), z)
display(fig)=#
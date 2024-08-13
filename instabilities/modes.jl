using Printf
using Statistics
include("diff-ops.jl")

function amplitude(u, v, w, b, f, aspect)
    
    return mean(abs.(u).^2) + mean(abs.(v).^2) + mean(abs.(w/aspect).^2) + mean(abs.(aspect/f * b).^2)

end

function least_stable_mode(Ri, k, l; H = 50, f = 1e-4, νᵥ = 1e-3, νₕ = 1e+1, N² = 1e-4, N = 100)
    
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
    u = 𝐯[0N+1:1N, index]
    v = 𝐯[1N+1:2N, index]
    w = 𝐯[2N+1:3N, index]
    b = 𝐯[3N+1:4N, index]
    p = 𝐯[4N+1:5N, index]

    A = amplitude(u, v, w, b, f, ((k^2+l^2)^0.5 * H)/2π)

    return u/A, v/A, w/A, b/A

end
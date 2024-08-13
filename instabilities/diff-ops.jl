using LinearAlgebra

function finite_diff_row(deriv_ord::Int; Δx::Real = 1, error_ord::Int = 2, offset::Int = 0)
    
    # Note that offset is taken from the LEFT-MOST gridpoint, not the centre
    L = deriv_ord
    h = error_ord
    n = offset
    # Finding fₙ⁽ᴸ⁾ with error 𝒪(Δxʰ)
    m = h + L   # Number of gridpoints to use (x₀ to xₘ₋₁)
    k_ind = -n : m-1-n
    Ã = [k^l / factorial(l) for k = k_ind, l = 0:m-1]   # Taylor series matrix
    #=row = OffsetMatrix(Ã^(-1), 0:m-1, k_ind)[L, :]      # Select correct row of Ã⁻¹
    # Offset so that 0 is the index of the place we want to evalute the derivative=#
    row = (Ã^(-1))[L+1, :]

    return round.(row / Δx^L, digits = 4)

end

function diff_mat(deriv_ord::Int, N::Int, Δx::Real; error_ord::Int = 2)

    if error_ord % 2 == 1
        println("ERROR: please use an error of even order. Increasing order of accuracy by 1")
        error_ord += 1
        # Central differences don't work unless error order is even
    end

    m = deriv_ord + error_ord
    # Minimum number of gridpoints required to calculate a derivative to this accuracy
    centred_offset = Int(floor((m-1)/2))
    D = zeros((N, N))
    for i = 1 : N
        distance_from_edge = minimum([i-1, N-i])
        if distance_from_edge > centred_offset
            row_short = finite_diff_row(deriv_ord, Δx = Δx, error_ord = error_ord, offset = centred_offset)
            row_short = (m % 2 == 0 ? row_short[1:m-1] : row_short) # Remove final redundant 0 if even number of grindpoints
            D[i, i-centred_offset : i+centred_offset] = row_short
        elseif distance_from_edge == i-1
            D[i, 1 : m] = finite_diff_row(deriv_ord, Δx = Δx, error_ord = error_ord, offset = i-1)
        else
            D[i, N-m+1 : N] = finite_diff_row(deriv_ord, Δx = Δx, error_ord = error_ord, offset = i-N+m-1)
        end
    end

    return D

end

function diff_mat(deriv_ord::Int; error_ord = 2)

    return diff_mat(deriv_ord, 10, 1.0, error_ord = error_ord)

end

function evals_with_constraints(A::Matrix, B::Matrix, C::Matrix, removed::Vector{Int64})
    
    # Solving the differential eigenvalue problem
    #           σA𝐯 = B𝐯
    # subject to contraints C𝐯 = 0.
    # removed is the set of indices of A and B where the
    # exact differential condition is relaxed

    N = size(A)[1]
    r = length(removed)
    k = N - r
    kept = Int.(zeros(k))
    j = 0
    for i = 1 : k
        j += 1
        while j in removed
            j += 1
        end
        kept[i] = j
    end
    Aₖₖ = A[kept, kept]
    Aₖᵣ = A[kept, removed]
    Bₖₖ = B[kept, kept]
    Bₖᵣ = B[kept, removed]
    Cᵣᵣ = C[:, removed]
    [@info Cᵣᵣ[i, :] for i in 1:size(Cᵣᵣ)[1]]
    Cᵣₖ = C[:, kept]

    G = - (Cᵣᵣ^(-1)) * Cᵣₖ
    Ã = Aₖₖ + Aₖᵣ * G
    B̃ = Bₖₖ + Bₖᵣ * G
    σ, 𝐯ₖ = eigen(B̃, Ã)
    𝐯ᵣ = G * 𝐯ₖ

    𝐯 = Complex.(zeros(N, k))
    𝐯[kept, :] = 𝐯ₖ
    𝐯[removed, :] = 𝐯ᵣ

    return σ, 𝐯

end

function Id_mat(N)

    return [i == j ? 1 : 0 for i = 1 : N, j = 1 : N]
    
end
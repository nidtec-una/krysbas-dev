function augmented_gram_schmidt_arnoldi(
    A,
    v::AbstractVector,
    m::Int,
    appendV::AbstractMatrix,
)
    n = size(A, 1)
    k = size(appendV, 2)
    s = m + k
    T = eltype(v)

    H = zeros(T, s + 1, s)
    V = zeros(T, n, s)
    w = zeros(T, n)       # single work vector (no n×s W matrix)

    copyto!(view(V, :, 1), v)

    for j = 1:s
        if j <= m
            mul!(w, A, view(V, :, j))
        else
            mul!(w, A, view(appendV, :, k - (j - m - 1)))
        end

        for i = 1:j
            vi = view(V, :, i)
            H[i, j] = dot(w, vi)
            axpy!(-H[i, j], vi, w)
        end
        h = norm(w)
        H[j+1, j] = h

        if h == 0
            return H[1:(j+1), 1:j], V[:, 1:j], j
        end
        # Near-breakdown: augmented direction is nearly in the Krylov space.
        # Drop the last (nearly-dependent) direction rather than normalising
        # w/h ≈ w/0 and producing a catastrophically ill-conditioned Rs.
        if j > 1 && h < sqrt(eps(T)) * H[2, 1]
            return H[1:j, 1:(j-1)], V[:, 1:(j-1)], j - 1
        end
        if j < s
            view(V, :, j + 1) .= w ./ h
        end
    end

    return H, V, s   # V is exactly n×s — no copy needed
end

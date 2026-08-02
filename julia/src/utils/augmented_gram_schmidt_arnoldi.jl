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
        w_norm = norm(w)   # pre-orthogonalisation scale, for the relative
        # near-breakdown check below -- NOT H[2, 1], which is anchored to
        # the unit-normalised fresh Krylov direction and does not shrink
        # as the algorithm converges, while an augmentation direction's
        # own scale does (it is built from the un-normalised correction
        # size). Comparing a shrinking quantity against that fixed
        # reference causes false-positive "breakdown" detection once the
        # residual gets small, silently disabling augmentation for many
        # cycles even though nothing is actually near-singular.

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
        # Near-breakdown: this direction is nearly in the existing basis.
        # Drop it rather than normalising w/h ≈ w/0 and producing a
        # catastrophically ill-conditioned Rs. Relative to the direction's
        # own pre-orthogonalisation norm, not an unrelated external scale.
        if j > 1 && h < sqrt(eps(T)) * w_norm
            return H[1:j, 1:(j-1)], V[:, 1:(j-1)], j - 1
        end
        if j < s
            view(V, :, j + 1) .= w ./ h
        end
    end

    return H, V, s   # V is exactly n×s — no copy needed
end

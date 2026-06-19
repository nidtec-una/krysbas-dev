function modified_gram_schmidt_arnoldi(A, v::AbstractVector, m::Int)
    n = size(A, 1)
    T = eltype(v)
    H = zeros(T, m + 1, m)
    V = zeros(T, n, m)
    w = zeros(T, n)       # single work vector (no n×m W matrix)

    copyto!(view(V, :, 1), v)

    for j = 1:m
        mul!(w, A, view(V, :, j))
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
        if j < m
            view(V, :, j + 1) .= w ./ h
        end
    end

    return H, V, m   # V is exactly n×m — no copy needed
end

function augmented_gram_schmidt_arnoldi(
        A, v::AbstractVector, m::Int, appendV::AbstractMatrix)
    n = size(A, 1)
    k = size(appendV, 2)
    s = m + k
    T = eltype(v)

    H = zeros(T, s + 1, s)
    W = zeros(T, n, s)
    V = zeros(T, n, s + 1)
    V[:, 1] = v

    for j in 1:s
        if j <= m
            W[:, j] = A * V[:, j]
        else
            # Augmented columns applied in reverse order (mirrors MATLAB fliplr logic)
            W[:, j] = A * appendV[:, k - (j - m - 1)]
        end

        for i in 1:j
            H[i, j] = dot(W[:, j], V[:, i])
            W[:, j] -= H[i, j] * V[:, i]
        end
        H[j + 1, j] = norm(W[:, j])

        if H[j + 1, j] == 0
            return H[1:j+1, 1:j], V[:, 1:j], j
        else
            V[:, j + 1] = W[:, j] / H[j + 1, j]
        end
    end

    return H, V[:, 1:s], s
end

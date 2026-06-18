function modified_gram_schmidt_arnoldi(A, v::AbstractVector, m::Int)
    n = size(A, 1)
    T = eltype(v)
    H = zeros(T, m + 1, m)
    W = zeros(T, n, m)
    V = zeros(T, n, m + 1)
    V[:, 1] = v

    for j in 1:m
        W[:, j] = A * V[:, j]
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

    return H, V[:, 1:m], m
end

function plane_rotations(H::AbstractMatrix, beta::Real)
    m = size(H, 2)

    g = zeros(eltype(H), m + 1)
    g[1] = beta

    for j in 1:m
        denom = hypot(H[j + 1, j], H[j, j])  # avoids sqrt(a²+b²) underflow when a,b are subnormal
        s = H[j + 1, j] / denom
        c = H[j, j] / denom

        P = Matrix{eltype(H)}(I, m + 1, m + 1)
        P[j, j] = c
        P[j + 1, j + 1] = c
        P[j, j + 1] = s
        P[j + 1, j] = -s

        H = P * H
        g = P * g
    end

    return H, g
end

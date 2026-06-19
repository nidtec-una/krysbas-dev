function plane_rotations(H::AbstractMatrix, beta::Real)
    m = size(H, 2)
    T = eltype(H)
    R = copy(H)
    g = zeros(T, m + 1)
    g[1] = beta

    for j in 1:m
        d = hypot(R[j, j], R[j + 1, j])
        c, s = d > 0 ? (R[j, j] / d, R[j + 1, j] / d) : (one(T), zero(T))

        # Apply Givens rotation in-place to the two affected rows (no P matrix)
        for k in j:m
            r_jk       =  c * R[j, k] + s * R[j + 1, k]
            R[j + 1, k] = -s * R[j, k] + c * R[j + 1, k]
            R[j, k]     = r_jk
        end
        g_j      =  c * g[j] + s * g[j + 1]
        g[j + 1] = -s * g[j] + c * g[j + 1]
        g[j]     = g_j
    end

    return R, g
end

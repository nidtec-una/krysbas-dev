function harmonic_ritz_vectors(F::AbstractMatrix, G::AbstractMatrix, k::Int,
        V::AbstractMatrix)
    # Solve the small dense generalised eigenvalue problem F*y = λ*G*y.
    # F and G are s×s (s = m+d, typically < 50), so eigen() is cheaper and
    # simpler than calling an iterative eigensolver like Arpack.
    vals, vecs = eigen(F, G)

    # Sort by ascending magnitude and take the k smallest.
    order = sortperm(abs.(vals))
    E = vecs[:, order[1:k]]

    # Lift the small eigenvectors back to the full space: yᵢ = V * eᵢ
    dy0 = V * E

    # If any eigenvector is complex, split into real and imaginary parts
    # (both are valid approximate eigenvectors of the real problem).
    if !isreal(dy0)
        dy = Matrix{Float64}(undef, size(dy0, 1), 0)
        for col in eachcol(dy0)
            if !isreal(col) && norm(real(col)) > 0
                dy = hcat(dy, real(col))
                if size(dy, 2) < k
                    dy = hcat(dy, abs.(imag(col)))
                end
            else
                dy = hcat(dy, real(col))
            end
            size(dy, 2) >= k && break
        end
        return dy
    end

    return real(dy0)
end

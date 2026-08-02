function harmonic_ritz_vectors(
    F::AbstractMatrix,
    G::AbstractMatrix,
    k::Int,
    V::AbstractMatrix,
)
    # Guard against a numerically indefinite G.
    #
    # Morgan's algorithm assumes G = R'*R is exactly SPD (it is, in exact
    # arithmetic, since R is the upper-triangular QR factor). After many
    # restart cycles, previously-added complex-pair augmentation vectors
    # (see below) can leave the augmented Krylov basis rank deficient, so
    # rounding breaks that assumption in practice and R develops a
    # near-zero diagonal entry, making G indefinite. Rather than let
    # eigen() error out on a non-positive-definite G, skip eigenvector
    # augmentation for this cycle: an empty dy makes the next call to
    # augmented_gram_schmidt_arnoldi fall back to a plain restarted GMRES
    # cycle, which restores a fresh, well-conditioned Krylov basis and
    # lets the eigenvector augmentation resume safely on a later cycle.
    #
    # Note this is NOT an ill-conditioning (cond) check: G = R'*R squares
    # the condition number of R by construction (this is exactly Morgan's
    # own G = R'*R shortcut, [1] eq. 16), so cond(G) is routinely as high
    # as 1e12-1e14 on perfectly healthy, converging cycles. A
    # relative-conditioning threshold was tried and triggered on nearly
    # every cycle as a false positive; only genuine loss of
    # positive-definiteness (isposdef failure) reliably distinguishes the
    # actual failure mode observed (eigen erroring on sherman5, d=5, in
    # the MATLAB/Octave port).
    if !isposdef(G)
        return zeros(eltype(V), size(V, 1), 0)
    end

    # Solve the small dense generalised eigenvalue problem F*y = λ*G*y.
    # F and G are s×s (s = m+d, typically < 50), so eigen() is cheaper and
    # simpler than calling an iterative eigensolver like Arpack.
    vals, vecs = eigen(F, G)

    # Sort ascending and take the k LARGEST magnitude.
    # Large |λ| in F*y = λ*G*y corresponds to small eigenvalues of A (the
    # harmonic Ritz values approximate eigenvalues from above), so we want LM.
    order = sortperm(abs.(vals))
    E = vecs[:, order[(end-k+1):end]]

    # Lift the small eigenvectors back to the full space: yᵢ = V * eᵢ
    dy0 = V * E

    # If any eigenvector is complex, split into real and imaginary parts
    # (both are valid approximate eigenvectors of the real problem).
    #
    # The inner append of the imaginary part is gated on the INPUT index ij
    # (matching MATLAB's harmonic_ritz_vectors.m exactly), not on the
    # OUTPUT column count so far: for k == 1, checking the output count
    # would read `1 < 1` right after the real part is appended and always
    # be false, silently dropping the imaginary part and returning only
    # half of a complex-conjugate pair -- the same class of bug found and
    # fixed in gmres_dr.m/.jl for exactly this k == 1, complex-pair case.
    if !isreal(dy0)
        dy = Matrix{Float64}(undef, size(dy0, 1), 0)
        ij = 1
        while size(dy, 2) <= k && ij <= k
            col = dy0[:, ij]
            if !isreal(col) && norm(real(col)) > 0
                dy = hcat(dy, real(col))
                if ij <= k
                    dy = hcat(dy, abs.(imag(col)))
                    if ij < k
                        ij = dy0[:, ij] == conj(dy0[:, ij+1]) ? ij + 2 : ij + 1
                    end
                end
            else
                dy = hcat(dy, real(col))
                ij += 1
            end
        end
        return dy
    end

    return real(dy0)
end

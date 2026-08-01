"""
    gmres_dr(A, b; m=0, k=-1, tol=1e-6, maxit=0, x_initial=[])

Restarted GMRES with deflated (thick) restarting (GMRES-DR(*m*, *k*)).

Unlike GMRES-E(*m*, *d*), which augments a fresh *m*-step Krylov subspace with
*d* extra eigenvectors (total dimension `m + d`), GMRES-DR(*m*, *k*) keeps the
subspace dimension fixed at *m* per cycle: the first `keep` basis vectors are
recycled harmonic Ritz vectors (`keep >= k`, see below) and the remaining
`m - keep` vectors come from standard Arnoldi.

The harmonic Ritz pairs are computed via the dense generalized eigenproblem
`F*y = λ*G*y` (`eigen(F, G)`), since these matrices are small (`m×m`). When
the `k` smallest-magnitude harmonic Ritz values include one half of a
complex-conjugate pair, both the real and imaginary parts of that eigenvector
are recycled instead of just one (equivalently `keep = k + 1`) — taking only
the real part would not actually lie in the associated 2-D real invariant
subspace and corrupts the thick restart's Arnoldi-consistency invariant.

The thick restart itself is a pure orthogonal basis rotation (Morgan 2002,
Wu & Simon 1999): the recycled vectors `Pk` and the current residual direction
are combined into an orthonormal `(keep+1)`-column basis `Pkp1` via thin QR,
and `H`, `V` are updated as `H <- Pkp1'*H*Pk`, `V <- V*Pkp1`. No matrix-vector
products with `A` are needed for the recycle step.

# Arguments
- `A`: square coefficient matrix (sparse or dense, `n×n`)
- `b::AbstractVector`: right-hand side vector of length `n`
- `m::Int=0`: subspace dimension per restart cycle; defaults to `min(n, 10)`.
  Setting `m == n` dispatches to full unrestarted GMRES.
- `k::Int=-1`: number of harmonic Ritz vectors to recycle at each restart;
  defaults to `min(m - 1, 3)`. Must satisfy `0 < k < m`. Setting `k == 0`
  dispatches to standard restarted GMRES(*m*).
- `tol::Real=1e-6`: relative residual tolerance for convergence
- `maxit::Int=0`: maximum number of restart cycles; defaults to `min(n, 10)`
- `x_initial::AbstractVector=[]`: initial guess; defaults to the zero vector

# Returns
- `x`: approximate solution vector
- `flag::Bool`: `true` if `relresvec[end] < tol` within `maxit` restarts
- `relresvec::Vector`: relative residual norm after each restart cycle
- `kdvec::Vector`: Krylov subspace dimension used at each cycle
- `time::Float64`: elapsed wall-clock time in seconds

# References
Morgan, R. B. (2002). GMRES with deflated restarting. *SIAM Journal on
Scientific Computing*, 24(1), 20-37.
[doi:10.1137/S1064827599364659](https://doi.org/10.1137/S1064827599364659)

Wu, K., & Simon, H. (2000). Thick-restart Lanczos method for large symmetric
eigenvalue problems. *SIAM Journal on Matrix Analysis and Applications*,
22(2), 602-616.
"""
function gmres_dr(
    A,
    b::AbstractVector;
    m::Int = 0,
    k::Int = -1,
    tol::Real = 1e-6,
    maxit::Int = 0,
    x_initial::AbstractVector = Float64[],
)

    # Sanity checks
    if ndims(A) != 2 || size(A, 1) == 0
        throw(ArgumentError("Matrix A cannot be empty."))
    end
    if size(A, 1) != size(A, 2)
        throw(ArgumentError("Matrix A must be square."))
    end
    n = size(A, 1)

    if isempty(b)
        throw(ArgumentError("Vector b cannot be empty."))
    end
    if length(b) != n
        throw(ArgumentError("Dimension mismatch between matrix A and vector b."))
    end

    if m == 0
        m = min(n, 10)
    end

    if isempty(x_initial)
        x_initial = zeros(eltype(b), n)
    end
    if length(x_initial) != n
        throw(
            ArgumentError(
                "Dimension mismatch between matrix A and initial guess x_initial.",
            ),
        )
    end

    # Dispatch: m == n → full unrestarted GMRES
    if m == n
        t0 = time()
        x, stats = Krylov.gmres(A, b; memory = n)
        elapsed = time() - t0
        res0 = norm(b - A * x_initial)
        resf = norm(b - A * x)
        relresvec = [1.0, resf / res0]
        kdvec = fill(n, 2)
        return x, stats.solved, relresvec, kdvec, elapsed
    end

    if m > n
        throw(ArgumentError("m must satisfy: 1 <= m <= n."))
    end

    # Set default k (-1 means unspecified)
    k = k == -1 ? min(m - 1, 3) : k

    # Dispatch: k == 0 explicitly → restarted GMRES(m)
    if k == 0
        t0 = time()
        x, stats = Krylov.gmres(A, b; memory = m)
        elapsed = time() - t0
        res0 = norm(b - A * x_initial)
        resf = norm(b - A * x)
        relresvec = [1.0, resf / res0]
        kdvec = fill(m, 2)
        return x, stats.solved, relresvec, kdvec, elapsed
    end

    if k >= m
        throw(ArgumentError("k must satisfy: 0 < k < m."))
    end

    eps_val = eps(Float64)
    if tol < eps_val
        @warn "Tolerance is too small; changed to eps."
        tol = eps_val
    elseif tol >= 1
        @warn "Tolerance is too large; changed to 1 - eps."
        tol = 1 - eps_val
    end

    if maxit == 0
        maxit = min(n, 10)
    end

    # --- GMRES-DR algorithm ---
    T = eltype(b)
    t0 = time()

    x = copy(x_initial)
    r0 = b - A * x
    beta = norm(r0)

    relresvec = [1.0]
    kdvec = Int[]

    V = zeros(T, n, m + 1)
    V[:, 1] = r0 / beta
    vr = [beta]
    H = zeros(T, 0, 0)
    keep = 0
    flag = false

    for _ = 1:maxit

        # Carry over the recycled (keep+1)-by-keep Hessenberg block's QR.
        q_qr = zeros(T, keep + 1, 0)
        r_qr = zeros(T, 0, 0)
        if keep > 0
            q_qr, r_qr = qrupdate_gs(H, q_qr, r_qr)
        end

        H_cycle = zeros(T, m + 1, m)
        H_cycle[1:size(H, 1), 1:size(H, 2)] .= H
        H = H_cycle

        # Arnoldi extension: add m - keep new basis vectors.
        converged = false
        for j = (keep+1):m
            w = A * V[:, j]
            for i = 1:j
                H[i, j] = dot(view(V, :, i), w)
                w .-= H[i, j] .* view(V, :, i)
            end
            H[j+1, j] = norm(w)
            V[:, j+1] = w / H[j+1, j]

            vr_c = vcat(vr, zeros(T, j + 1 - length(vr)))
            q_qr, r_qr = qrupdate_gs(view(H, 1:(j+1), 1:j), q_qr, r_qr)

            d = r_qr \ (q_qr' * vr_c)
            rc = vr_c - H[1:(j+1), 1:j] * d
            res = norm(rc)

            if res < tol * beta
                x = x + V[:, 1:j] * d
                push!(relresvec, res / beta)
                push!(kdvec, j)
                flag = true
                converged = true
                break
            end
        end
        converged && return x, flag, relresvec, kdvec, time() - t0

        # Finalise this cycle's solution/residual from the completed QR
        # factorisation. Needed even when the loop above did not execute at
        # all (keep == m: the harmonic Ritz selection consumed the entire
        # per-cycle budget, leaving no room for a fresh Arnoldi step).
        vr_c = vcat(vr, zeros(T, m + 1 - length(vr)))
        d = r_qr \ (q_qr' * vr_c)
        rc = vr_c - H[1:(m+1), 1:m] * d
        res = norm(rc)

        x = x + V[:, 1:m] * d
        push!(relresvec, res / beta)
        push!(kdvec, m)

        if relresvec[end] < tol
            flag = true
            return x, flag, relresvec, kdvec, time() - t0
        end

        # --- Harmonic Ritz: dense generalised eigenproblem F*y = λ*G*y ---
        g_mat = r_qr' * r_qr
        f_mat = H[1:m, 1:m]'

        vals, vecs = eigen(f_mat, g_mat)
        # F*y = λ*G*y is the reciprocal formulation (Morgan 2000, eq. 2.10):
        # λ = 1/θ̃, so the smallest-magnitude harmonic Ritz values θ̃ (the
        # ones we want to deflate) correspond to the LARGEST-magnitude λ
        # here. Walk descending |λ|, i.e. ascending order reversed.
        order = reverse(sortperm(abs.(vals)))

        # Collect k real recycling directions. A complex eigenvalue
        # contributes both the real and imaginary parts of its eigenvector
        # (the 2-D real invariant subspace of its conjugate pair); its
        # partner is then marked consumed so it is not processed twice.
        consumed = falses(m)
        cols = Vector{T}[]
        for idx in order
            length(cols) >= k && break
            consumed[idx] && continue
            lambda = vals[idx]
            v = vecs[:, idx]
            consumed[idx] = true
            if isreal(v)
                push!(cols, real(v))
            else
                push!(cols, real(v))
                push!(cols, imag(v))
                for j in order
                    if !consumed[j] &&
                       isapprox(vals[j], conj(lambda); atol = 1e-8 * max(abs(lambda), 1))
                        consumed[j] = true
                        break
                    end
                end
            end
        end
        keep = length(cols)
        pk = Matrix(qr(reduce(hcat, cols)).Q)[:, 1:keep]

        # Thick restart: orthogonal basis rotation (see docstring).
        pkp1_raw = hcat(vcat(pk, zeros(T, 1, keep)), rc)
        pkp1 = Matrix(qr(pkp1_raw).Q)[:, 1:(keep+1)]
        pkp1[1:m, 1:keep] = pk   # restore exact sign of recycled columns

        H = pkp1' * H[1:(m+1), 1:m] * pk
        V[:, 1:(keep+1)] = V[:, 1:(m+1)] * pkp1
        vr = pkp1' * rc
    end

    return x, flag, relresvec, kdvec, time() - t0
end

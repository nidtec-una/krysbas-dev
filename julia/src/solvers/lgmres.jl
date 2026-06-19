"""
    lgmres(A, b; m=0, l=-1, tol=1e-6, maxit=0, x_initial=[])

Restarted GMRES augmented with error approximation vectors (LGMRES(*m*, *l*)).

At each restart, up to `l` error approximation vectors from prior cycles are
appended to the Krylov subspace, preserving information that standard restarted
GMRES discards. This typically yields faster convergence than GMRES(*m*+*l*)
at the same per-cycle cost.

# Arguments
- `A`: square coefficient matrix (sparse or dense, `n×n`)
- `b::AbstractVector`: right-hand side vector of length `n`
- `m::Int=0`: base Krylov restart dimension; defaults to `min(n, 10)`. Setting
  `m == n` dispatches to full unrestarted GMRES.
- `l::Int=-1`: maximum number of stored error vectors to append; defaults to `3`.
  Setting `l == 0` dispatches to standard restarted GMRES(*m*).
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
Baker, A. H., Jessup, E. R., & Manteuffel, T. (2005). A technique for
accelerating the convergence of restarted GMRES. *SIAM Journal on Matrix
Analysis and Applications*, 26(4), 962–984.
[doi:10.1137/S0895479803422014](https://doi.org/10.1137/S0895479803422014)
"""
function lgmres(A, b::AbstractVector;
                m::Int=0,
                l::Int=-1,
                tol::Real=1e-6,
                maxit::Int=0,
                x_initial::AbstractVector=Float64[])

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
    if m > n
        throw(ArgumentError("m must satisfy: 1 <= m <= n."))
    end

    if isempty(x_initial)
        x_initial = zeros(eltype(b), n)
    end
    if length(x_initial) != n
        throw(ArgumentError("Dimension mismatch between matrix A and initial guess x_initial."))
    end

    # Dispatch: m == n → full unrestarted GMRES
    if m == n
        t0 = time()
        x, stats = Krylov.gmres(A, b; memory=n)
        elapsed = time() - t0
        r0 = b - A * x_initial
        res0 = norm(r0)
        resf = norm(b - A * x)
        relresvec = [1.0, resf / res0]
        kdvec = fill(n, 2)
        return x, stats.solved, relresvec, kdvec, elapsed
    end

    # Dispatch: l == 0 explicitly → restarted GMRES(m)
    if l == 0
        t0 = time()
        x, stats = Krylov.gmres(A, b; memory=m)
        elapsed = time() - t0
        r0 = b - A * x_initial
        res0 = norm(r0)
        resf = norm(b - A * x)
        relresvec = [1.0, resf / res0]
        kdvec = fill(m, 2)
        return x, stats.solved, relresvec, kdvec, elapsed
    end

    # Set default l (-1 means unspecified)
    l = l == -1 ? 3 : l

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

    T = eltype(b)
    t0 = time()

    r_buf = b - A * x_initial   # residual buffer, reused each restart
    res0 = norm(r_buf)
    relresvec = [1.0]

    # Sliding window of error approximation vectors (oldest → newest: col 1 → col l)
    z_mat = zeros(T, n, l)
    z_current = zeros(T, n)   # reused each restart to avoid per-cycle allocation

    # Iteration 1: GMRES(m+l) cycle with online Givens rotations and convergence
    # detection within the cycle — mirrors MATLAB's gmres(A,b,m+l,tol,1,...,x0).
    s_max = m + l
    H1 = zeros(T, s_max + 1, s_max)
    V1 = zeros(T, n, s_max + 1)
    c1 = zeros(T, s_max)    # Givens cosines
    sv1 = zeros(T, s_max)   # Givens sines
    gv1 = zeros(T, s_max + 1)
    gv1[1] = res0
    beta = res0
    view(V1, :, 1) .= r_buf ./ beta
    s1 = s_max  # actual subspace dim for this cycle

    for j in 1:s_max
        wj = view(V1, :, j + 1)   # use next V1 column as work buffer
        mul!(wj, A, view(V1, :, j))
        for i in 1:j
            vi = view(V1, :, i)
            H1[i, j] = dot(wj, vi)
            axpy!(-H1[i, j], vi, wj)
        end
        h_norm = norm(wj)
        H1[j + 1, j] = h_norm
        for i in 1:j - 1  # apply previous rotations to column j
            tmp           =  c1[i] * H1[i, j] + sv1[i] * H1[i + 1, j]
            H1[i + 1, j]  = -sv1[i] * H1[i, j] + c1[i]  * H1[i + 1, j]
            H1[i, j]      = tmp
        end
        denom = hypot(H1[j, j], H1[j + 1, j])
        if denom > 0
            c1[j]  = H1[j, j]     / denom
            sv1[j] = H1[j + 1, j] / denom
        else
            c1[j] = one(T); sv1[j] = zero(T)
        end
        H1[j, j] = denom; H1[j + 1, j] = zero(T)
        gv1[j + 1] = -sv1[j] * gv1[j]
        gv1[j]     =   c1[j] * gv1[j]
        s1 = j
        if h_norm > 0 && j < s_max
            wj ./= h_norm   # normalize in place; wj IS V1[:,j+1]
        end
        if abs(gv1[j + 1]) / res0 < tol || h_norm == 0
            break
        end
    end

    Rs1 = UpperTriangular(H1[1:s1, 1:s1])
    minimizer = Rs1 \ gv1[1:s1]
    x_cur = copy(x_initial)
    mul!(x_cur, view(V1, :, 1:s1), minimizer, 1.0, 1.0)   # x_cur = x_initial + V1*min

    push!(relresvec, abs(gv1[s1 + 1]) / res0)
    mul!(view(z_mat, :, 1), view(V1, :, 1:s1), minimizer)  # z1 = x - x_initial

    if relresvec[end] < tol
        kdvec = fill(m + l, length(relresvec))
        return x_cur, true, relresvec, kdvec, time() - t0
    end

    flag = false

    # Main loop: restart = 2, 3, ...
    for restart in 2:maxit
        copyto!(r_buf, b)
        mul!(r_buf, A, x_cur, -1.0, 1.0)   # r_buf = b - A*x_cur
        beta = norm(r_buf)
        r_buf ./= beta                       # r_buf = v1 (normalized)

        n_z = min(restart - 1, l)           # stored error vectors to use
        m_eff = m + l - n_z                  # effective Krylov dimension

        z_slice = view(z_mat, :, 1:n_z)
        H, V, s = augmented_gram_schmidt_arnoldi(A, r_buf, m_eff, z_slice)
        HUpTri, g = plane_rotations(H, beta)

        Rs = HUpTri[1:s, 1:s]
        minimizer = Rs \ g[1:s]

        # z_current = V[:,1:m_eff_actual]*min[1:m_eff_actual]
        #           + z_slice[:,reversed]*min[m_eff_actual+1:s]
        # Avoids copy(V) by splitting the matvec across the two sub-matrices.
        m_eff_actual = min(m_eff, s)
        n_aug = s - m_eff_actual
        mul!(z_current, view(V, :, 1:m_eff_actual), view(minimizer, 1:m_eff_actual))
        for i in 1:n_aug
            axpy!(minimizer[m_eff_actual + i],
                  view(z_slice, :, n_z - i + 1), z_current)
        end

        x_cur .+= z_current
        push!(relresvec, abs(g[s + 1]) / res0)

        if relresvec[end] < tol
            flag = true
            kdvec = fill(m + l, length(relresvec))
            return x_cur, flag, relresvec, kdvec, time() - t0
        end

        # Update sliding window of error vectors
        if restart <= l
            copyto!(view(z_mat, :, restart), z_current)
        else
            for i in 1:l - 1
                copyto!(view(z_mat, :, i), view(z_mat, :, i + 1))
            end
            copyto!(view(z_mat, :, l), z_current)
        end
    end

    kdvec = fill(m + l, length(relresvec))
    return x_cur, flag, relresvec, kdvec, time() - t0
end

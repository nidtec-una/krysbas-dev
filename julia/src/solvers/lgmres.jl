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

    r0 = b - A * x_initial
    res0 = norm(r0)
    relresvec = [1.0]

    # Sliding window of error approximation vectors (oldest → newest: col 1 → col l)
    z_mat = zeros(T, n, l)

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
    V1[:, 1] = r0 / beta
    s1 = s_max  # actual subspace dim for this cycle

    for j in 1:s_max
        w = A * V1[:, j]
        for i in 1:j
            H1[i, j] = dot(w, V1[:, i])
            w -= H1[i, j] * V1[:, i]
        end
        h_norm = norm(w)
        H1[j + 1, j] = h_norm
        for i in 1:j - 1  # apply previous rotations to column j
            tmp         =  c1[i] * H1[i, j] + sv1[i] * H1[i + 1, j]
            H1[i + 1, j] = -sv1[i] * H1[i, j] + c1[i]  * H1[i + 1, j]
            H1[i, j] = tmp
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
            V1[:, j + 1] = w / h_norm
        end
        if abs(gv1[j + 1]) / res0 < tol || h_norm == 0
            break
        end
    end

    Rs1 = UpperTriangular(H1[1:s1, 1:s1])
    minimizer = Rs1 \ gv1[1:s1]
    x = x_initial + V1[:, 1:s1] * minimizer

    push!(relresvec, abs(gv1[s1 + 1]) / res0)
    z_mat[:, 1] = x - x_initial

    if relresvec[end] < tol
        kdvec = fill(m + l, length(relresvec))
        return x, true, relresvec, kdvec, time() - t0
    end

    x_cur = x
    flag = false

    # Main loop: restart = 2, 3, ...
    for restart in 2:maxit
        r = b - A * x_cur
        beta = norm(r)
        v1 = r / beta

        n_z = min(restart - 1, l)     # stored error vectors to use
        m_eff = m + l - n_z            # effective Krylov dimension

        z_slice = z_mat[:, 1:n_z]     # error vectors (oldest col 1 … newest col n_z)
        H, V, s = augmented_gram_schmidt_arnoldi(A, v1, m_eff, z_slice)
        HUpTri, g = plane_rotations(H, beta)

        Rs = HUpTri[1:s, 1:s]
        minimizer = Rs \ g[1:s]

        # W: Krylov columns from V, augmented columns from actual error vectors
        # (newest error vector first, matching fliplr(zMat[:,1:n_z]) in MATLAB)
        m_eff_actual = min(m_eff, s)
        W = copy(V)  # already n × s
        n_aug = s - m_eff_actual
        if n_aug > 0
            W[:, m_eff_actual + 1:s] = z_slice[:, n_z:-1:n_z - n_aug + 1]
        end

        z_current = W * minimizer
        xm = x_cur + z_current

        push!(relresvec, abs(g[s + 1]) / res0)

        if relresvec[end] < tol
            flag = true
            kdvec = fill(m + l, length(relresvec))
            return xm, flag, relresvec, kdvec, time() - t0
        end

        x_cur = xm

        # Update sliding window of error vectors
        if restart <= l
            z_mat[:, restart] = z_current
        else
            z_mat[:, 1:l - 1] = z_mat[:, 2:l]
            z_mat[:, l] = z_current
        end
    end

    kdvec = fill(m + l, length(relresvec))
    return x_cur, flag, relresvec, kdvec, time() - t0
end

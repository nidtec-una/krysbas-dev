function gmres_e(A, b::AbstractVector;
                 m::Int=0,
                 d::Int=-1,
                 tol::Real=1e-6,
                 maxit::Int=0,
                 x_initial::AbstractVector=Float64[],
                 eigstol::Real=1e-6)

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

    # Dispatch: d == 0 explicitly → restarted GMRES(m)
    if d == 0
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

    # Set default d (-1 means unspecified)
    d = d == -1 ? min(m, 3) : d
    if d > m
        throw(ArgumentError("d cannot be larger than m"))
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

    # GMRES-E algorithm
    t0 = time()
    flag = false
    r0 = b - A * x_initial
    res0 = norm(r0)
    relresvec = [1.0]
    beta = res0
    v1 = r0 / beta

    # Iteration 1: plain MGS Arnoldi with s = m + d augmented dimension
    H, V, s = modified_gram_schmidt_arnoldi(A, v1, m + d)
    HUpTri, g = plane_rotations(H, beta)

    Rs = HUpTri[1:s, 1:s]
    minimizer = Rs \ g[1:s]
    xm = x_initial + V * minimizer

    push!(relresvec, abs(g[s + 1]) / res0)

    if relresvec[end] < tol
        flag = true
        return xm, flag, relresvec, fill(s, length(relresvec)), time() - t0
    end

    # Harmonic Ritz vectors for first augmentation (Morgan 1995, p.1161 step 5)
    Fold = H[1:s, 1:s]'
    G = Rs' * Rs
    dy = harmonic_ritz_vectors(Fold, G, d, V)

    # Main GMRES-E loop (outer iterations 2, 3, ...)
    x = xm
    for restart in 2:maxit
        r = b - A * xm
        beta = norm(r)
        v1 = r / beta

        H, V, s = augmented_gram_schmidt_arnoldi(A, v1, m, dy[:, 1:d])
        HUpTri, g = plane_rotations(H, beta)

        Rs = HUpTri[1:s, 1:s]
        minimizer = Rs \ g[1:s]

        # Replace augmented columns with harmonic Ritz vectors (step 4 in [1])
        n_aug = max(0, s - m)
        V_use = copy(V)
        if n_aug > 0
            V_use[:, m + 1:s] = dy[:, 1:n_aug]
        end

        x = xm + V_use * minimizer
        push!(relresvec, abs(g[s + 1]) / res0)

        if relresvec[end] < tol
            flag = true
            return x, flag, relresvec, fill(s, length(relresvec)), time() - t0
        elseif restart < maxit
            W = V_use[:, 1:s]
            Fold = W' * (A' * W)
            G = Rs' * Rs
            dy = harmonic_ritz_vectors(Fold, G, d, V_use)
        end

        xm = x
    end

    return x, flag, relresvec, fill(s, length(relresvec)), time() - t0
end

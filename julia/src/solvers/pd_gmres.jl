function pd_gmres(A, b::AbstractVector;
                  m_initial::Int=0,
                  m_min_max::Union{Vector{Int},Nothing}=nothing,
                  m_step::Int=1,
                  tol::Real=1e-6,
                  maxit::Int=0,
                  x_initial::AbstractVector=Float64[],
                  alpha_pd::Vector{Float64}=[-3.0, 5.0])

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

    if isempty(x_initial)
        x_initial = zeros(eltype(b), n)
    end
    if length(x_initial) != n
        throw(ArgumentError("Dimension mismatch between matrix A and initial guess x_initial."))
    end

    # m_initial == 0 means "not given" → unrestarted; also m_initial == n is unrestarted
    restarted = m_initial != 0 && m_initial != n

    if restarted && (m_initial < 1 || m_initial > n)
        throw(ArgumentError("m_initial must satisfy: 1 <= m_initial <= n."))
    end

    # m_min_max defaults
    m_min, m_max = 1, n
    if !isnothing(m_min_max)
        if !restarted
            @warn "m_min_max was given but will not be used."
        else
            m_min, m_max = m_min_max[1], m_min_max[2]
            if m_min < 1 || m_max > n || m_max <= m_min
                throw(ArgumentError("m_min_max must satisfy: 1 <= m_min_max[1] < m_min_max[2] <= n."))
            end
            if m_min > m_initial || m_max < m_initial
                throw(ArgumentError("m_min_max must satisfy: m_min_max[1] <= m_initial <= m_min_max[2]."))
            end
        end
    end

    if m_step < 1 || m_step >= n
        throw(ArgumentError("m_step must satisfy: 0 < m_step < n."))
    end

    eps_val = eps(Float64)
    if tol < eps_val
        @warn "Tolerance is too small; changed to eps."
        tol = eps_val
    elseif tol >= 1
        @warn "Tolerance is too large; changed to 1 - eps."
        tol = 1 - eps_val
    end

    alpha_p, alpha_d = alpha_pd[1], alpha_pd[2]

    # Unrestarted: solve A*(x-x_initial) = b-A*x_initial, then recover x.
    if !restarted
        if maxit == 0
            maxit = min(n, 10)
        end
        t0 = time()
        r0_unr = b - A * x_initial
        res0 = norm(r0_unr)
        x_shift, stats = Krylov.gmres(A, r0_unr; memory=n, itmax=n)
        x = x_initial + x_shift
        elapsed = time() - t0
        resf = norm(b - A * x)
        relresvec = [1.0, resf / res0]
        kdvec = [NaN]
        return x, resf / res0 < tol, relresvec, kdvec, elapsed
    end

    # Restarted PD-GMRES
    if maxit == 0
        maxit = min(ceil(Int, n / m_initial), 10)
    end

    t0 = time()
    flag = false
    m_cur = m_initial
    x_cur = copy(x_initial)

    r0 = b - A * x_cur
    res0 = norm(r0)
    res_abs = [res0]  # absolute residuals for pd_rule
    relresvec = [1.0]
    kdvec = [m_initial]

    for restart in 1:typemax(Int)
        # Update m via PD rule (not on first iteration)
        if restart > 1
            m_cur, m_initial = pd_rule(m_cur, n, m_initial, m_min, m_max, m_step,
                                       res_abs, restart, alpha_p, alpha_d)
        end
        push!(kdvec, m_cur)

        r = b - A * x_cur
        beta = norm(r)
        v1 = r / beta

        H, V, m_used = modified_gram_schmidt_arnoldi(A, v1, m_cur)
        HUpTri, g = plane_rotations(H, beta)

        Rs = HUpTri[1:m_used, 1:m_used]
        minimizer = Rs \ g[1:m_used]
        xm = x_cur + V * minimizer

        push!(res_abs, abs(g[m_used + 1]))
        push!(relresvec, res_abs[end] / res0)

        if relresvec[end] < tol || length(relresvec) >= maxit
            flag = relresvec[end] < tol
            return xm, flag, relresvec, kdvec, time() - t0
        end

        x_cur = xm
        m_cur = m_used  # may be less than requested if happy breakdown occurred
    end

    return x_cur, flag, relresvec, kdvec, time() - t0
end

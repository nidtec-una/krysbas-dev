"""
    pd_gmres(A, b; m_initial=0, m_min_max=nothing, m_step=1,
             tol=1e-6, maxit=0, x_initial=[], alpha_pd=[-3.0, 5.0])

Restarted GMRES with a Proportional-Derivative controller (PD-GMRES(*m*)).

The restart parameter *m* is adjusted automatically each cycle by a PD control
law driven by the recent residual history, avoiding the need to hand-tune a
fixed restart value.

# Arguments
- `A`: square coefficient matrix (sparse or dense, `n×n`)
- `b::AbstractVector`: right-hand side vector of length `n`
- `m_initial::Int=0`: initial restart dimension. `0` or `n` selects unrestarted
  GMRES; any value in `1:n-1` enables the PD-adaptive restarted mode.
- `m_min_max::Union{Vector{Int},Nothing}=nothing`: two-element vector
  `[m_min, m_max]` bounding the adaptive restart range; ignored when unrestarted.
  Defaults to `[1, n]`.
- `m_step::Int=1`: increment/decrement step size for the PD rule
- `tol::Real=1e-6`: relative residual tolerance for convergence
- `maxit::Int=0`: maximum number of restart cycles; defaults to
  `ceil(n / m_initial)` in restarted mode, or `min(n, 10)` unrestarted
- `x_initial::AbstractVector=[]`: initial guess; defaults to the zero vector
- `alpha_pd::Vector{Float64}=[-3.0, 5.0]`: PD gain vector `[αₚ, α_d]`

# Returns
- `x`: approximate solution vector
- `flag::Bool`: `true` if `relresvec[end] < tol` within `maxit` restarts
- `relresvec::Vector`: relative residual norm after each restart cycle
- `kdvec::Vector`: restart dimension *m* used at each cycle
- `time::Float64`: elapsed wall-clock time in seconds

# References
Núñez, R. C., Schaerer, C. E., & Bhaya, A. (2018). A proportional-derivative
control strategy for restarting the GMRES(*m*) algorithm. *Journal of
Computational and Applied Mathematics*, 337, 209–224.
[doi:10.1016/j.cam.2018.01.009](https://doi.org/10.1016/j.cam.2018.01.009)
"""
function pd_gmres(
    A,
    b::AbstractVector;
    m_initial::Int = 0,
    m_min_max::Union{Vector{Int},Nothing} = nothing,
    m_step::Int = 1,
    tol::Real = 1e-6,
    maxit::Int = 0,
    x_initial::AbstractVector = Float64[],
    alpha_pd::Vector{Float64} = [-3.0, 5.0],
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
                throw(
                    ArgumentError(
                        "m_min_max must satisfy: 1 <= m_min_max[1] < m_min_max[2] <= n.",
                    ),
                )
            end
            if m_min > m_initial || m_max < m_initial
                throw(
                    ArgumentError(
                        "m_min_max must satisfy: m_min_max[1] <= m_initial <= m_min_max[2].",
                    ),
                )
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
        x_shift, stats = Krylov.gmres(A, r0_unr; memory = n, itmax = n)
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

    r_buf = b - A * x_cur   # residual buffer, reused each restart
    res0 = norm(r_buf)
    res_abs = [res0]  # absolute residuals for pd_rule
    relresvec = [1.0]
    kdvec = [m_initial]

    for restart = 1:typemax(Int)
        # Update m via PD rule (not on first iteration)
        if restart > 1
            m_cur, m_initial = pd_rule(
                m_cur,
                n,
                m_initial,
                m_min,
                m_max,
                m_step,
                res_abs,
                restart,
                alpha_p,
                alpha_d,
            )
        end
        push!(kdvec, m_cur)

        copyto!(r_buf, b)
        mul!(r_buf, A, x_cur, -1.0, 1.0)   # r_buf = b - A*x_cur
        beta = norm(r_buf)
        r_buf ./= beta                       # r_buf = v1 (normalized)

        H, V, m_used = modified_gram_schmidt_arnoldi(A, r_buf, m_cur)
        HUpTri, g = plane_rotations(H, beta)

        Rs = HUpTri[1:m_used, 1:m_used]
        minimizer = Rs \ g[1:m_used]
        mul!(x_cur, V, minimizer, 1.0, 1.0)   # x_cur += V * minimizer

        push!(res_abs, abs(g[m_used+1]))
        push!(relresvec, res_abs[end] / res0)

        if relresvec[end] < tol || length(relresvec) >= maxit
            flag = relresvec[end] < tol
            return x_cur, flag, relresvec, kdvec, time() - t0
        end

        m_cur = m_used  # may be less than requested if happy breakdown occurred
    end

    return x_cur, flag, relresvec, kdvec, time() - t0
end

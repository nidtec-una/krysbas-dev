"""
    a_slgmres_e(A, b; m_initial=0, m_min_max=nothing, m_step=1, l=-1, d=-1,
                epsilon_threshold=0.01, alpha_pd=[2.0, 0.8], tol=1e-6, maxit=0,
                x_initial=[], eigstol=1e-6)

A-SLGMRES-E(*mⱼ*, *l*, *d*): extends [`slgmres_e`](@ref) with an adaptive
restart parameter.

Whenever a cycle triggers the same convergence-slowdown signal that switches
`slgmres_e` from LGMRES-style to GMRES-E-style augmentation, this variant
ALSO grows the restart parameter *m* via the proportional-derivative law of
`pd_rule` (the same law used by [`pd_gmres`](@ref)), for as long as slowdown
persists. *m* only changes on stagnating cycles; it stays fixed during
LGMRES-style cycles. This is Algorithm 1 of [1].

See `slgmres_e` for the switching rule itself (eq. 33-35 of [1]) and why it
can be read off the Givens-rotated residual at no extra cost.

Deviation from the reference implementation: `pd_rule`'s "warm-up" gating
(it only applies the derivative term once at least 3 cycles of residual
history exist, and only the proportional term with 2) is driven here by the
complete, correctly-indexed cycle history. Cabral's own reference script
(`Adaptive_PD_lgmres_e.m`) has a residual-history update for its first cycle
commented out, which leaves its own cycle counter one cycle behind for the
rest of the run. Both this port and the reference converge correctly; they
simply grow *m* on a slightly different schedule.

# Arguments
- `A`: square coefficient matrix (sparse or dense, `n×n`)
- `b::AbstractVector`: right-hand side vector of length `n`
- `m_initial::Int=0`: initial restart parameter; defaults to `min(n, 10)`.
  Setting `m_initial == n` dispatches to full unrestarted GMRES.
- `m_min_max::Union{Vector{Int},Nothing}=nothing`: two-element vector
  `[m_min, m_max]` bounding *m*'s adaptive range. Defaults to `[1, n]`.
- `m_step::Int=1`: step size for growing the internal floor `pd_rule` falls
  back on when the PD law proposes `m < m_min`.
- `l::Int=-1`: number of error approximation vectors to append during
  LGMRES-style cycles; defaults to `3`. Must satisfy `l > 0`.
- `d::Int=-1`: number of harmonic Ritz vectors to append during GMRES-E
  -style cycles; defaults to `min(m_initial, 3)`. Must satisfy `d > 0`.
- `epsilon_threshold::Real=0.01`: slowdown threshold. A cycle is classified
  as stagnating when `‖r^(j)‖ / ‖r^(j-1)‖ >= 1 - epsilon_threshold`.
- `alpha_pd::Vector{Float64}=[2.0, 0.8]`: proportional and derivative
  coefficients for the restart-parameter growth law (`pd_rule`). These are
  the values reported in [1] -- note this differs from `pd_gmres`'s own
  default of `[-3.0, 5.0]`, which comes from a different paper ([2]) tuned
  for a different purpose (shrinking as well as growing *m* every cycle,
  rather than only growing it on detected stagnation).
- `tol::Real=1e-6`: relative residual tolerance for convergence
- `maxit::Int=0`: maximum number of restart cycles; defaults to `min(n, 10)`
- `x_initial::AbstractVector=[]`: initial guess; defaults to the zero vector
- `eigstol::Real=1e-6`: reserved for API compatibility (unused internally,
  since the underlying `harmonic_ritz_vectors` uses a dense eigensolver
  here, not an iterative one)

# Returns
- `x`: approximate solution vector
- `flag::Bool`: `true` if `relresvec[end] < tol` within `maxit` restarts
- `relresvec::Vector`: relative residual norm after each restart cycle
- `kdvec::Vector`: Krylov subspace dimension used at each cycle (`m` during
  cycle 1, then `m+l` or `m+d` thereafter -- with `m` itself possibly having
  grown -- depending on which augmentation was active)
- `time::Float64`: elapsed wall-clock time in seconds

# References
[1] Cabral, J. C., Schaerer, C. E., & Bhaya, A. (2020). Improving GMRES(*m*)
using an adaptive switching controller. *Numerical Linear Algebra with
Applications*, 27(5), e2305.
[doi:10.1002/nla.2305](https://doi.org/10.1002/nla.2305)

[2] Núñez, R. C., Schaerer, C. E., & Bhaya, A. (2018). A
proportional-derivative control strategy for restarting the GMRES(*m*)
algorithm. *Journal of Computational and Applied Mathematics*, 337, 209–224.
[doi:10.1016/j.cam.2018.01.009](https://doi.org/10.1016/j.cam.2018.01.009)
"""
function a_slgmres_e(
    A,
    b::AbstractVector;
    m_initial::Int = 0,
    m_min_max::Union{Vector{Int},Nothing} = nothing,
    m_step::Int = 1,
    l::Int = -1,
    d::Int = -1,
    epsilon_threshold::Real = 0.01,
    alpha_pd::Vector{Float64} = [2.0, 0.8],
    tol::Real = 1e-6,
    maxit::Int = 0,
    x_initial::AbstractVector = Float64[],
    eigstol::Real = 1e-6,
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

    if m_initial == 0
        m_initial = min(n, 10)
    end

    # Dispatch: m_initial == n → full unrestarted GMRES
    if m_initial == n
        t0 = time()
        x, stats = Krylov.gmres(A, b; memory = n)
        elapsed = time() - t0
        res0 = norm(b - A * x_initial)
        resf = norm(b - A * x)
        relresvec = [1.0, resf / res0]
        kdvec = fill(n, 2)
        return x, stats.solved, relresvec, kdvec, elapsed
    end

    if m_initial < 1 || m_initial > n
        throw(ArgumentError("m_initial must satisfy: 1 <= m_initial <= n."))
    end

    m_min, m_max = 1, n
    if !isnothing(m_min_max)
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

    if m_step < 1 || m_step > n - 1
        throw(ArgumentError("m_step must satisfy: 0 < m_step < n."))
    end

    l = l == -1 ? 3 : l
    d = d == -1 ? min(m_initial, 3) : d

    if l <= 0
        throw(ArgumentError("l must satisfy: l > 0."))
    end
    if d <= 0
        throw(ArgumentError("d must satisfy: d > 0."))
    end

    if epsilon_threshold <= 0 || epsilon_threshold >= 1
        throw(ArgumentError("epsilon_threshold must satisfy: 0 < epsilon_threshold < 1."))
    end

    alpha_p, alpha_d = alpha_pd[1], alpha_pd[2]

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

    # --- A-SLGMRES-E algorithm ---
    T = eltype(b)
    norm_y = 1 - epsilon_threshold

    x = copy(x_initial)
    r0 = b - A * x
    res1 = norm(r0)

    relresvec = [1.0]
    kdvec = Int[]

    # Current restart parameter and pd_rule's own "initial" bookkeeping
    # variable (see pd_rule.jl / pd_gmres.jl: m_current tracks the growth
    # floor pd_rule falls back on when the PD law proposes m < m_min).
    m = m_initial
    m_current = m_initial

    # Sliding window of LGMRES-style error approximation vectors, newest
    # vector last. n_z tracks how many columns are actually populated.
    z_mat = zeros(T, n, l)
    n_z = 0

    t0 = time()

    # -------------------------------------------------------------------
    # Cycle 1: plain GMRES(m). Neither augmentation strategy, nor the
    # restart-parameter growth law (which needs cycle history), applies
    # yet.
    # -------------------------------------------------------------------
    v1 = r0 / res1
    H, V, s = modified_gram_schmidt_arnoldi(A, v1, m)
    h_up_tri, g = plane_rotations(H, res1)

    rs = h_up_tri[1:s, 1:s]
    gs = g[1:s]
    minimizer = rs \ gs
    z_cycle = V * minimizer
    x = x + z_cycle
    push!(relresvec, abs(g[s+1]) / res1)
    push!(kdvec, s)

    if relresvec[end] < tol
        return x, true, relresvec, kdvec, time() - t0
    end

    z_mat[:, 1] = z_cycle
    n_z = 1

    stagnating = relresvec[end] / relresvec[end-1] >= norm_y
    dy = zeros(T, n, 0)
    if stagnating
        # Cycle 1's basis is a plain (non-augmented) Arnoldi basis, so the
        # cheap H'-based formula for fold is exact here -- see slgmres_e.jl
        # for why this shortcut is only valid on this first cycle.
        fold = H[1:s, 1:s]'
        g_mat = rs' * rs
        dy = harmonic_ritz_vectors(fold, g_mat, d, V)
    end

    # -------------------------------------------------------------------
    # Main loop: cycles 2, 3, ...
    # -------------------------------------------------------------------
    flag = false
    while !flag && length(relresvec) - 1 < maxit

        # ------------------------------------------------------------------
        # Control block: on a stagnating cycle, grow m via the PD law
        # (pd_rule, Algorithm 1 of [2]) before building this cycle's
        # subspace. m is left untouched on non-stagnating (LGMRES-style)
        # cycles -- growth only happens in direct response to detected
        # slowdown, matching [1]'s Adaptive_PD_lgmres_e.m reference.
        # ------------------------------------------------------------------
        if stagnating
            m, m_current = pd_rule(
                m,
                n,
                m_current,
                m_min,
                m_max,
                m_step,
                relresvec,
                length(relresvec),
                alpha_p,
                alpha_d,
            )
        end

        r = b - A * x
        beta = norm(r)
        v1 = r / beta

        local s_cyc
        if !stagnating
            # --- LGMRES-style cycle ---
            # See slgmres_e.jl for why the reverse() placement here is the
            # opposite of the GMRES-E-style branch below.
            l_use = min(n_z, l)
            processing_order = reverse(z_mat[:, 1:l_use], dims = 2)
            H, V, s_cyc = augmented_gram_schmidt_arnoldi(A, v1, m, z_mat[:, 1:l_use])
            h_up_tri, g = plane_rotations(H, beta)
            rs = h_up_tri[1:s_cyc, 1:s_cyc]
            gs = g[1:s_cyc]
            minimizer = rs \ gs
            V[:, (m+1):s_cyc] = processing_order[:, 1:(s_cyc-m)]
        else
            # --- GMRES-E-style cycle (using the just-grown m) ---
            #
            # dy IS sliced to d columns here, matching slgmres_e.jl: this
            # family of algorithms is built on Cabral's reference
            # implementations, whose GMRES-E-style branch (e.g.
            # Adaptive_PD_lgmres_e.m) hard-codes s = m + d and only ever
            # reads dy[:, 1:d] -- functionally identical to this explicit
            # slice, even though harmonic_ritz_vectors can return more
            # columns when a harmonic Ritz value is complex. See
            # slgmres_e.jl's GMRES-E-style branch for the full rationale,
            # and gmres_e.jl (which has no such reference to match) for
            # why that solver instead uses dy in full.
            processing_order = dy[:, 1:d]
            H, V, s_cyc = augmented_gram_schmidt_arnoldi(
                A,
                v1,
                m,
                reverse(processing_order, dims = 2),
            )
            h_up_tri, g = plane_rotations(H, beta)
            rs = h_up_tri[1:s_cyc, 1:s_cyc]
            gs = g[1:s_cyc]
            minimizer = rs \ gs
            V[:, (m+1):s_cyc] = processing_order[:, 1:(s_cyc-m)]
        end
        s = s_cyc

        aux = V * minimizer
        x = x + aux
        push!(relresvec, abs(g[s+1]) / res1)
        push!(kdvec, s)

        if relresvec[end] < tol
            flag = true
            break
        end

        # ------------------------------------------------------------------
        # Decide the augmentation strategy for the NEXT cycle -- see
        # slgmres_e.jl for the full rationale, which carries over unchanged
        # here (the growth of m does not affect this decision).
        # ------------------------------------------------------------------
        ratio = relresvec[end] / relresvec[end-1]

        if ratio >= norm_y
            stagnating = true
            W = V[:, 1:s]
            fold = W' * A' * W
            g_mat = rs' * rs
            dy = harmonic_ritz_vectors(fold, g_mat, d, V)
        else
            stagnating = false
            if n_z < l
                n_z += 1
                z_mat[:, n_z] = aux
            else
                z_mat[:, 1:(l-1)] = z_mat[:, 2:l]
                z_mat[:, l] = aux
            end
        end
    end

    return x, flag, relresvec, kdvec, time() - t0
end

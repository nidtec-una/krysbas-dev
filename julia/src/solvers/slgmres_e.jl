"""
    slgmres_e(A, b; m=0, l=-1, d=-1, epsilon_threshold=0.01, tol=1e-6, maxit=0, x_initial=[], eigstol=1e-6)

Restarted GMRES that switches, cycle by cycle, between two augmentation
strategies based on an observed convergence-slowdown signal (SLGMRES-E(*m*,
*l*, *d*)).

By default each cycle behaves like LGMRES(*m*, *l*): the restart subspace is
augmented with up to *l* error approximation vectors from prior cycles. When
slowdown is detected -- the relative residual ratio between two consecutive
cycles exceeds `1 - epsilon_threshold` -- it switches to GMRES-E(*m*, *d*)
for one or more cycles: the subspace is instead augmented with *d* harmonic
Ritz vectors computed from the cycle that triggered the detection. It
switches back to LGMRES-style augmentation as soon as the ratio improves
again.

The slowdown signal, `epsilon` in `‖r_m^(j)‖ / ‖r_m^(j-1)‖ = 1 - epsilon`, is
read directly off the Givens-rotated residual already computed for the
least-squares solve every cycle -- no extra cost. Under this switching rule
the residual norm is provably non-increasing cycle to cycle (Theorem 2 of
[1]).

# Arguments
- `A`: square coefficient matrix (sparse or dense, `n×n`)
- `b::AbstractVector`: right-hand side vector of length `n`
- `m::Int=0`: restart parameter, fixed across all cycles; defaults to
  `min(n, 10)`. Setting `m == n` dispatches to full unrestarted GMRES.
- `l::Int=-1`: number of error approximation vectors to append during
  LGMRES-style cycles; defaults to `3`. Must satisfy `l > 0`.
- `d::Int=-1`: number of harmonic Ritz vectors to append during GMRES-E
  -style cycles; defaults to `min(m, 3)`. Must satisfy `d > 0`.
- `epsilon_threshold::Real=0.01`: slowdown threshold. A cycle is classified
  as stagnating when `‖r^(j)‖ / ‖r^(j-1)‖ >= 1 - epsilon_threshold`. Default
  matches the numerical experiments of [1].
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
  cycle 1, `m+l` or `m+d` thereafter depending on which augmentation was
  active)
- `time::Float64`: elapsed wall-clock time in seconds

# References
Cabral, J. C., Schaerer, C. E., & Bhaya, A. (2020). Improving GMRES(m) using
an adaptive switching controller. *Numerical Linear Algebra with
Applications*, 27(5), e2305.
[doi:10.1002/nla.2305](https://doi.org/10.1002/nla.2305)
"""
function slgmres_e(
    A,
    b::AbstractVector;
    m::Int = 0,
    l::Int = -1,
    d::Int = -1,
    epsilon_threshold::Real = 0.01,
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

    l = l == -1 ? 3 : l
    d = d == -1 ? min(m, 3) : d

    if l <= 0
        throw(ArgumentError("l must satisfy: l > 0."))
    end
    if d <= 0
        throw(ArgumentError("d must satisfy: d > 0."))
    end

    if epsilon_threshold <= 0 || epsilon_threshold >= 1
        throw(ArgumentError("epsilon_threshold must satisfy: 0 < epsilon_threshold < 1."))
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

    # --- SLGMRES-E algorithm ---
    T = eltype(b)
    norm_y = 1 - epsilon_threshold

    x = copy(x_initial)
    r0 = b - A * x
    res1 = norm(r0)

    relresvec = [1.0]
    kdvec = Int[]

    # Sliding window of LGMRES-style error approximation vectors, newest
    # vector last. n_z tracks how many columns are actually populated
    # (independent of how many cycles have elapsed, since GMRES-E-style
    # cycles do not touch this window).
    z_mat = zeros(T, n, l)
    n_z = 0

    t0 = time()

    # -------------------------------------------------------------------
    # Cycle 1: plain GMRES(m). Neither augmentation strategy has any
    # history to draw on yet.
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

    # The correction from this cycle is always kept as the first error
    # approximation vector, regardless of what the next cycle turns out
    # to be.
    z_mat[:, 1] = z_cycle
    n_z = 1

    stagnating = relresvec[end] / relresvec[end-1] >= norm_y
    dy = zeros(T, n, 0)
    if stagnating
        # Cycle 1's basis is a plain (non-augmented) Arnoldi basis, so the
        # cheap H'-based formula for fold is exact here (equivalent to,
        # but cheaper than, W'*A'*W). Later cycles cannot use this
        # shortcut because V's augmented columns get overwritten with raw
        # (non-Arnoldi-consistent) vectors below -- see the main loop.
        fold = H[1:s, 1:s]'
        g_mat = rs' * rs
        dy = harmonic_ritz_vectors(fold, g_mat, d, V)
    end

    # -------------------------------------------------------------------
    # Main loop: cycles 2, 3, ...
    # -------------------------------------------------------------------
    flag = false
    while !flag && length(relresvec) - 1 < maxit

        r = b - A * x
        beta = norm(r)
        v1 = r / beta

        local s_cyc
        if !stagnating
            # --- LGMRES-style cycle ---
            # Note the reverse() placement here is the opposite of the
            # GMRES-E-style branch below: matching lgmres.jl's own
            # convention (not gmres_e.jl's), the newest error vector goes
            # into the first augmentation slot V[:, m+1]. processing_order
            # is that same slot order (newest first); augmented_gram_schmidt
            # _arnoldi processes it in this order and, unlike the MATLAB
            # utility, may drop trailing (nearly-dependent) columns via its
            # own near-breakdown guard, so the V overwrite below must not
            # assume all l_use columns survived.
            l_use = min(n_z, l)
            processing_order = reverse(z_mat[:, 1:l_use], dims = 2)
            H, V, s_cyc = augmented_gram_schmidt_arnoldi(A, v1, m, z_mat[:, 1:l_use])
            h_up_tri, g = plane_rotations(H, beta)
            rs = h_up_tri[1:s_cyc, 1:s_cyc]
            gs = g[1:s_cyc]
            minimizer = rs \ gs
            V[:, (m+1):s_cyc] = processing_order[:, 1:(s_cyc-m)]
        else
            # --- GMRES-E-style cycle ---
            processing_order = dy[:, 1:d]
            H, V, s_cyc =
                augmented_gram_schmidt_arnoldi(A, v1, m, reverse(processing_order, dims = 2))
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
        # Decide the augmentation strategy for the NEXT cycle from the
        # ratio just observed. This mirrors [1] eq. (33)-(35): the ratio is
        # read off the residual we already computed above, at no extra
        # cost. The decision is independent of which strategy produced THIS
        # cycle's result -- a recovering GMRES-E cycle's correction is just
        # as valid an error-approximation vector as an LGMRES cycle's, and
        # a newly-stagnating LGMRES cycle's basis is just as valid a source
        # of harmonic Ritz vectors as a GMRES-E cycle's.
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

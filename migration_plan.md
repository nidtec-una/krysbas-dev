# Migration Plan: KrySBAS MATLAB → Julia

## Status

**Steps 0–10 (GMRES-E, LGMRES, PD-GMRES) are complete** — merged via PR #81
(2026-06-19).

**GMRES-DR, SLGMRES-E, and A-SLGMRES-E have since been added to the MATLAB
side and ported to Julia as well** (Steps 11–13 below), on the `switching`
branch. The original plan below is kept as-written for historical reference;
new steps follow the same format.

## Dependency analysis

**Standard GMRES** is well-implemented in Julia. **LGMRES, GMRES-E, and PD-GMRES are not in
any Julia package** — all three must be ported. However, the fallback dispatch cases (`m == n`
→ unrestarted GMRES, `d == 0` / `l == 0` → standard restarted GMRES) can delegate to
`Krylov.jl` instead of reimplementing them.

| Dependency | Replaces | Notes |
|---|---|---|
| `Krylov.jl` | MATLAB built-in `gmres` | Fallback dispatch in all three solvers; the LGMRES first-cycle call |
| `LinearAlgebra` (stdlib) | `eigs`, `norm`, `\`, etc. | `eigen(F, G)` replaces `eigs(F, G, k, 'LM')` — valid because F, G are small dense matrices (size `m+d`) |
| `SparseArrays` (stdlib) | MATLAB sparse | Transparent to algorithm logic |
| `MAT.jl` | `load('data/X.mat')` | Test data only — `.mat` files stay as-is |
| `Test` (stdlib) | MOxUnit | `@test`, `@test_throws`, `@testset` |

**No `Arpack.jl` needed.** The generalized eigenvalue problem in `harmonic_ritz_vectors`
operates on matrices of size `s×s` where `s = m+d` (typically < 30). `LinearAlgebra.eigen(F,
G)` computes all eigenvalues of a small dense matrix, which is cleaner and faster than calling
Arpack for a handful of eigenpairs.

**`plane_rotations` can be ported directly.** Julia's `LinearAlgebra.givens` exists but the
naive loop (building rotation matrix P and multiplying) is fine at this subspace size, and a
direct port minimises translation risk.

---

## Package structure

Original plan (Steps 0–10):

```
KrySBAS.jl/
├── Project.toml
├── src/
│   ├── KrySBAS.jl            # module, includes, exports
│   ├── solvers/
│   │   ├── gmres_e.jl
│   │   ├── lgmres.jl
│   │   └── pd_gmres.jl
│   └── utils/
│       ├── plane_rotations.jl
│       ├── modified_gram_schmidt_arnoldi.jl
│       ├── augmented_gram_schmidt_arnoldi.jl
│       ├── harmonic_ritz_vectors.jl
│       └── pd_rule.jl
└── test/
    ├── runtests.jl
    ├── test_plane_rotations.jl
    ├── test_modified_gram_schmidt_arnoldi.jl
    ├── test_augmented_gram_schmidt_arnoldi.jl
    ├── test_harmonic_ritz_vectors.jl
    ├── test_pd_rule.jl
    ├── test_gmres_e.jl
    ├── test_lgmres.jl
    ├── test_pd_gmres.jl
    └── test_poisson.jl
```

Current structure, after Steps 11–13 (GMRES-DR, SLGMRES-E, A-SLGMRES-E):

```
KrySBAS.jl/
├── Project.toml
├── src/
│   ├── KrySBAS.jl
│   ├── solvers/
│   │   ├── gmres_e.jl
│   │   ├── lgmres.jl
│   │   ├── pd_gmres.jl
│   │   ├── gmres_dr.jl
│   │   ├── slgmres_e.jl
│   │   └── a_slgmres_e.jl
│   └── utils/
│       ├── plane_rotations.jl
│       ├── modified_gram_schmidt_arnoldi.jl
│       ├── augmented_gram_schmidt_arnoldi.jl
│       ├── harmonic_ritz_vectors.jl
│       ├── pd_rule.jl
│       └── qrupdate_gs.jl        # GMRES-DR's QR update step
└── test/
    ├── runtests.jl
    ├── test_plane_rotations.jl
    ├── test_modified_gram_schmidt_arnoldi.jl
    ├── test_augmented_gram_schmidt_arnoldi.jl
    ├── test_harmonic_ritz_vectors.jl
    ├── test_pd_rule.jl
    ├── test_gmres_e.jl
    ├── test_lgmres.jl
    ├── test_pd_gmres.jl
    ├── test_gmres_dr.jl
    ├── test_qrupdate_gs.jl
    ├── test_slgmres_e.jl
    ├── test_a_slgmres_e.jl
    └── test_poisson.jl
```

---

## Translation notes that apply across the whole codebase

These differences appear in every file:

- **`nargin` → default arguments.** Every `if (nargin < k) || isempty(x); x = default; end`
  becomes a Julia keyword argument with a default: `function foo(A, b; m=min(size(A,1), 10),
  ...)`. Use `nothing` when a default must be computed from other args (e.g. `maxit` depends on
  both `n` and `m`).
- **`error()` / `warning()` → same names in Julia**, different string syntax. Use `error("...")`
  or `@warn "..."`.
- **MATLAB `eps` → Julia `eps(Float64)`.**
- **`tic/toc` → `time()`** captured as `t0 = time(); ...; elapsed = time() - t0`, or use the
  `@elapsed` macro.
- **`size(A)` returns a tuple in Julia** — `n, _ = size(A)` works; `[n, ~] = size(A)` does not
  exist.
- **`isempty` exists in Julia** and behaves the same.
- **`zeros(n, 1)` → `zeros(n)`.** Julia vectors are 1D (`Vector{Float64}`), not column
  matrices. Prefer `zeros(n)` and `Vector` throughout.
- **`\` for backslash solve** works the same.
- **Broadcasting:** MATLAB's `a .* b` → Julia's `a .* b` (same syntax, but Julia's
  dot-broadcasting is more general).
- **`fliplr(v)` → `reverse(v, dims=2)`** for matrices, or `reverse(v)` for vectors.

---

## Step-by-step plan

Each step ends with a green test run before moving on.

---

### Step 0 — Package scaffold and CI skeleton

1. Create the package with `] generate KrySBAS` in Julia.
2. Write `Project.toml` with `[deps]` for `Krylov` and `MAT`, and `[compat]` bounds.
3. Write `src/KrySBAS.jl`: empty module with `include` stubs and placeholder exports.
4. Write `test/runtests.jl` that `include`s each test file.
5. Write a minimal GitHub Actions workflow (`.github/workflows/julia_tests.yml`) using
   `julia-actions/setup-julia` + `julia -e 'using Pkg; Pkg.test("KrySBAS")'`.
6. Write a `julia_style.yml` using `JuliaFormatter.jl` to enforce formatting (replaces
   MISS_HIT).

**Test checkpoint:** `Pkg.test()` runs with zero test files — no errors, no output. CI passes
on an empty module.

---

### Step 1 — `plane_rotations`

Port `plane_rotations.m` directly. The signature becomes:

```julia
function plane_rotations(H::Matrix, beta::Real)::Tuple{Matrix, Vector}
```

Key differences: `[~, m] = size(H)` → `m = size(H, 2)`. The loop body is otherwise identical.

**Tests:** Port `plane_rotations` tests verifying:
- A known 2×1 Hessenberg example gives expected upper-triangular output.
- `g[1] == beta` before rotations.
- Output is type-stable (both outputs are `Float64`).

---

### Step 2 — `modified_gram_schmidt_arnoldi`

Direct port. Signature:

```julia
function modified_gram_schmidt_arnoldi(A, v::Vector, m::Int)
    # returns (H, V, m_updated)
```

The early-exit on `H[j+1, j] == 0` (happy breakdown) stays identical.

**Tests:**
- Identity matrix + unit vector → H and V have expected structure.
- Happy breakdown triggers early return with truncated `m`.
- Output V has orthonormal columns (`V' * V ≈ I`).

---

### Step 3 — `augmented_gram_schmidt_arnoldi`

Direct port. Signature:

```julia
function augmented_gram_schmidt_arnoldi(A, v::Vector, m::Int, appendV::Matrix)
    # returns (H, V, s)
```

The logic around the augmented block (`j <= m` vs. `j > m`) is identical.

**Tests:**
- With `appendV` of zero columns (`k=0`) should match `modified_gram_schmidt_arnoldi` output.
- Augmented columns increase `s` correctly.
- `V` columns are orthonormal.

---

### Step 4 — `harmonic_ritz_vectors`

This is the only utility that needs a substantive change. Replace:

```matlab
[E2, D2] = eigs(F, G, k, 'LM', opts);
```

with:

```julia
vals, vecs = eigen(F, G)   # all eigenvalues of small dense problem
```

Then sort `abs.(vals)` ascending and take the first `k` eigenvectors. This replaces the
`'LM'` + sort-ascending pattern in the MATLAB code with a cleaner all-eigenvalues-then-select
approach (valid because `s` is always small).

The complex-conjugate splitting logic ports directly — `isreal(v)` → `!isreal(v)`, and
`real(v)`, `imag(v)` all exist in Julia.

**Tests:**
- Known small (3×3, 4×4) generalized eigenvalue problem with precomputed harmonic Ritz vectors.
- Complex eigenvector case: verify split produces two real vectors with correct norms.
- Output has `k` columns.

---

### Step 5 — `pd_rule`

Pure arithmetic, direct port. Signature:

```julia
function pd_rule(m, n, mInitial, mMin, mMax, mStep,
                 res::Vector, iter::Int, alphaP, alphaD)::Tuple{Int, Int}
```

Return a `(mj, mInitial)` tuple instead of the `[mj mInitial]` row vector.

**Tests:** Port existing numerical tests with known `res` vectors — verify output `m` matches
expected values for the `iter > 3`, `iter > 2`, and `iter <= 2` branches.

---

### Step 6 — `gmres_e`

Signature:

```julia
function gmres_e(A, b::Vector;
                 m::Int = min(size(A, 1), 10),
                 d::Int = -1,           # -1 signals "use default min(m, 3)"
                 tol::Float64 = 1e-6,
                 maxit::Int = min(size(A, 1), 10),
                 x0::Union{Vector, Nothing} = nothing,
                 eigstol::Float64 = 1e-6)
```

Fallback dispatch (replacing MATLAB's built-in `gmres`):

```julia
# m == n: unrestarted
x, stats = Krylov.gmres(A, b; atol=0.0, rtol=tol, itmax=n)

# m < n, d == 0: standard restarted GMRES
x, stats = Krylov.gmres(A, b; restart=m, atol=0.0, rtol=tol, itmax=maxit * m)
```

The rest of the algorithm (Arnoldi loop, plane rotations, harmonic Ritz update) calls the
utility functions from Steps 1–4.

**Tests:** Port all tests from `test_gmres_e.m`:
- Input validation (`@test_throws`).
- Fallback cases (identity matrix, `m == n`, `d == 0`).
- Embree 3×3 toy example.
- Sherman1, Sherman4 sparse matrices loaded via `MAT.load("data/sherman1.mat")`.

---

### Step 7 — `lgmres`

The LGMRES first-cycle call is currently `gmres(A, b, m+l, tol, 1, ...)` (one cycle of
GMRES(m+l)). In Julia, replace with:

```julia
x, stats = Krylov.gmres(A, b; restart=m + l, atol=0.0, rtol=tol,
                         itmax=m + l, x0=x0)
```

The rest (augmented Arnoldi loop, `zMat` history, error vector accumulation) is a direct port
of the while loop.

**Tests:** Port all tests from `test_lgmres.m` — same matrix suite as GMRES-E.

---

### Step 8 — `pd_gmres`

The unrestarted path:

```julia
x, stats = Krylov.gmres(A, b; atol=0.0, rtol=tol, itmax=maxit, x0=x0)
```

The restarted path uses `modified_gram_schmidt_arnoldi` + `plane_rotations` + `pd_rule` already
ported above — direct translation of the while loop.

**Tests:** Port all tests from `test_pd_gmres.m`.

---

### Step 9 — Integration test

Port `test_poisson.m`. This exercises all three solvers on a real PDE-derived system. Verify:
- All three solvers converge on the Poisson problem.
- `relresvec[end] < tol` for each solver.

---

### Step 10 — Final CI and formatter

- Update `julia_tests.yml` to run the full test suite including the integration test.
- Run `JuliaFormatter.format(".")` and commit formatted code.
- Add `codecov` upload step (Julia coverage via `Pkg.test` with `coverage=true`).
- Remove the MATLAB workflows, or keep them temporarily if a parallel maintenance period is
  needed.

---

### Step 11 — `qrupdate_gs` and `gmres_dr`

GMRES-DR(*m*, *k*) restarts with harmonic Ritz vectors recovered from a Schur
decomposition of the previous cycle's Hessenberg matrix, then extends the
Krylov basis with a QR-update step (`qrupdate_gs`) instead of rebuilding the
augmented Arnoldi basis from scratch each cycle. Both `qrupdate_gs.m` and
`gmres_dr.m` port directly — same `nargin`/`[~, m] = size(...)`/`fliplr`
translation patterns as Steps 1–8.

One correctness fix surfaced during porting and was applied to **both**
languages: `dy`/harmonic-Ritz-vector handling for a complex-conjugate
harmonic Ritz pair must keep the vector count in whatever's actually
returned, not implicitly assume `k`/`d` columns — silently truncating to
`d` discards the imaginary half of a legitimate conjugate pair. See
`gmres_dr.m`/`.jl` and `qrupdate_gs.m`/`.jl`.

**Julia-only addition:** `augmented_gram_schmidt_arnoldi.jl` gained a
relative near-breakdown check (`h < sqrt(eps(T)) * w_norm`, dropping the
offending column) beyond MATLAB's exact-zero-only check. This was needed to
fix a real, confirmed `lgmres.jl` divergence from `lgmres.m` (87 vs. 38
cycles on sherman1) caused by an unrelated redundant augmentation direction
going undetected. Porting the equivalent check to MATLAB was tried and
reverted — it caused a regression in `lgmres.m`'s own embree3 test (a false
near-breakdown positive on a tiny, exactly-representable 3×3 system) with no
compensating benefit, since MATLAB's `lgmres.m` was not exhibiting the bug
the check was meant to fix. The two languages are intentionally *not*
symmetric here.

**Tests:** Port all tests from `test_gmres_dr.m`/`test_qrupdate_gs.m`. The
embree3 toy example documents a known, expected stall (GMRES-DR(2,1) is the
only non-trivial (*m*, *k*) combination for a 3×3 system, and the complex
conjugate Ritz pair there forces `keep = k+1 = m`, leaving no room for a
fresh Arnoldi direction) rather than asserting convergence.

---

### Step 12 — `slgmres_e`

SLGMRES-E(*m*, *l*, *d*) (Cabral, Schaerer & Bhaya, 2020) switches, cycle by
cycle, between LGMRES-style and GMRES-E-style augmentation based on a
convergence-slowdown signal read directly off the Givens-rotated residual —
no extra cost. It's a genuinely new solver on the MATLAB side too (not a
pre-existing file), implemented directly from the published algorithm and
validated against Cabral's own reference scripts (`jcc_codigos_may_2023/`,
not committed to this repo) before porting to Julia.

Direct port of `slgmres_e.m`'s structure: cycle 1 is plain GMRES(*m*); the
main loop's two branches reuse `augmented_gram_schmidt_arnoldi` exactly as
`lgmres.jl` and `gmres_e.jl` already do, just switched per cycle. One
non-obvious point ported faithfully in both directions: the LGMRES-style
branch's augmentation columns go into `V[:, m+1:s]` in newest-first order
(matching `lgmres.jl`'s own convention), while the GMRES-E-style branch's go
in natural order (matching `gmres_e.jl`'s convention) — the two branches use
*opposite* `reverse`/`fliplr` placements at the call site for this reason.

**`dy` truncation policy — deliberately different from `gmres_e`:**
`gmres_e.m`/`.jl` pass the full harmonic-Ritz-vector matrix through (no
truncation to `d` columns), matching Morgan (1995) step 5 literally, since
there is no external reference implementation to match and an apparent
"truncation converges faster" result under benchmarking turned out to be
noise from a chaotic stagnation boundary rather than a real effect (see
`gmres_e.m`'s in-code comment for the full writeup). `slgmres_e.m`/`.jl`
keep the `dy(:, 1:d)` truncation instead, because — unlike `gmres_e` — this
solver *is* validated cycle-by-cycle against an authoritative reference
(Cabral's own scripts), and that reference truncates too (implicitly, via a
hard-coded `s = m+d` loop bound in `Adaptive_lgmres_e_switch.m` that only
ever reads `dy(:, 1:d)`). Match the ground truth you actually have.

`harmonic_ritz_vectors.m`/`.jl` also gained a guard against `G` losing
positive-definiteness before the generalized eigenproblem solve (a real,
pre-existing crash independent of the truncation question — e.g. `gmres_e`
on sherman5 with `d=5` errors with `eigs: matrix B is not positive
definite` even without any dy-truncation change). On failure it skips
augmentation for that cycle rather than erroring, self-healing on the next
plain restart. An initial attempt paired this with an `rcond`/`cond`-based
ill-conditioning threshold; that had to be dropped because `G = R'*R`
squares `R`'s condition number by construction (Morgan's own shortcut, [1]
eq. 16), so `rcond(G)` is routinely `1e-12`–`1e-14` on perfectly healthy
cycles — only a true loss of positive-definiteness (`chol`/`isposdef`
failure) reliably distinguishes the real failure mode.

**Tests:** Port all tests from `test_slgmres_e.m`. Its embree3 toy example
is a Julia-only structural stall: at this size, both augmentation branches
propose a direction that's redundant up to machine precision
(`h/w_norm ~ 2e-16`), which Julia's near-breakdown check (see Step 11)
correctly declines, permanently capping the subspace at `s = m`. MATLAB's
unguarded Arnoldi proceeds anyway and, numerically lucky on this exact
3×3 system, still lands on the right answer — the same divergence class
documented for GMRES-DR's embree3 test in Step 11.

---

### Step 13 — `a_slgmres_e`

Extends `slgmres_e` with a PD-adaptive restart parameter: on a stagnating
cycle, *m* grows via `pd_rule` (Algorithm 1 of [1]) before that cycle's
subspace is built, reusing `pd_rule.jl` and its existing `pd_gmres.jl`
calling convention (`(mj, m_initial) = pd_rule(m, n, m_initial, m_min,
m_max, m_step, res, iter, alpha_p, alpha_d)`) rather than introducing a new
one.

**Known, documented deviation from the reference** (both MATLAB and Julia):
`pd_rule`'s warm-up gating is driven by the complete, correctly-indexed
cycle history. Cabral's own reference script
(`Adaptive_PD_lgmres_e.m`) has a residual-history update for its first cycle
commented out, leaving its own cycle counter one cycle behind for the rest
of the run. Both this port and the reference converge correctly on the
sherman5 case from [1]'s own numerical experiments; they simply grow *m* on
a slightly different schedule (MATLAB port: 118 cycles, max *m*+*d* = 62;
Julia port: 131 cycles, max *m*+*d* = 60 — same ballpark, not exact, and not
expected to be).

**Tests:** Port all tests from `test_a_slgmres_e.m`. Same embree3 structural
stall as `slgmres_e.jl` (Step 12), for the same reason.

References:
[1] Cabral, J. C., Schaerer, C. E., & Bhaya, A. (2020). Improving GMRES(*m*)
using an adaptive switching controller. *Numerical Linear Algebra with
Applications*, 27(5), e2305.

---

## Installing Julia locally

```bash
# 1. Install juliaup (version manager)
brew install juliaup
juliaup add release      # installs latest stable Julia

# 2. VS Code extension: "Julia" (julialang.language-julia)

# 3. In the Julia REPL — add development dependencies
]  # enters package mode
add Krylov MAT JuliaFormatter
```

Running tests locally at any step:

```julia
] test KrySBAS             # from the package root in pkg mode
# or equivalently
using Pkg; Pkg.test()
```

Running a single test file during development:

```julia
include("test/test_gmres_e.jl")
```

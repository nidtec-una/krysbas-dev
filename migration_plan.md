# Migration Plan: KrySBAS MATLAB → Julia

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

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

KrySBAS is a MATLAB toolbox of adaptive iterative solvers for sparse linear systems (Ax = b), built on Krylov subspaces. The three solvers are:

- **GMRES-E(m, d)** — restarted GMRES augmented with `d` harmonic Ritz vectors (approximate eigenvectors of smallest eigenvalues)
- **LGMRES(m, l)** — restarted GMRES augmented with `l` error approximation vectors from prior cycles
- **PD-GMRES(m)** — restarted GMRES with a Proportional-Derivative controller that adapts the restart parameter `m` automatically

All solvers share the same output signature: `[x, flag, relresvec, kdvec, time]`.

## Commands

### Lint (MISS_HIT)
```bash
make lint                  # runs mh_style --fix && mh_metric --ci && mh_lint inside matlab/
cd matlab && mh_style --fix
cd matlab && mh_metric --ci
cd matlab && mh_lint
```

### Run tests with coverage (requires MATLAB and MOxUnit/MOcov installed)
```bash
make coverage      # runs MATLAB -r "run_tests; exit()" inside matlab/, opens matlab/coverage_html/
```

To run tests directly in MATLAB:
```matlab
% From matlab/:
run run_tests()

% Run a single test file:
addpath(genpath('src'));
moxunit_runtests('tests/test_gmres_e.m', '-verbose')
```

### Install test dependencies (for local MATLAB)
Clone MOxUnit and MOcov into the repo root (they are gitignored by the CI workflow):
```bash
git clone https://github.com/MOxUnit/MOxUnit.git --depth 1
git clone https://github.com/MOcov/MOcov.git --depth 1
```

## Code style rules (matlab/miss_hit.cfg)

- Line length: **80 characters**
- Function names: `snake_case` — regex `[a-zA-Z]+(_[a-zA-Z0-9]+)*`
- Script names: `snake_case` — regex `[a-zA-Z]+(_[a-zA-Z]+)*`
- Parameter names: `snake_case` — regex `[a-zA-Z]+(_[a-zA-Z]+)*`
- Max nesting depth (`cnest`): 6
- Max file length: 500 lines
- Max cyclomatic complexity (`cyc`): 15
- Max parameters: 5
- Copyright header required: `CC&MA - NIDTec - FP - UNA`
- The `data/` directory is excluded from linting

## Architecture

```
matlab/
  src/
    solvers/        # The three solver implementations
      gmres_e.m
      lgmres.m
      pd_gmres.m
    utils/          # Shared mathematical building blocks
      modified_gram_schmidt_arnoldi.m   # MGS Arnoldi (standard)
      augmented_gram_schmidt_arnoldi.m  # Arnoldi augmented with vectors (for GMRES-E/LGMRES)
      harmonic_ritz_vectors.m           # Eigenvalue problem for augmentation vectors
      pd_rule.m                         # PD control law for restart parameter update
      plane_rotations.m                 # Givens rotations → upper triangular H + g vector
  tests/            # MOxUnit test files, one per solver plus Poisson integration test
  run_tests.m       # Entry point for test runner with coverage
  miss_hit.cfg      # MISS_HIT style/metric config (project_root set here)
julia/
  src/              # Julia package (KrySBAS.jl) — see migration_plan.md
  test/
data/               # Shared .mat files with sparse test matrices (embree3, sherman1/4/5)
```

**Solver internals follow this pattern:**
1. Sanity checks on all inputs with explicit error messages
2. Special-case dispatch (e.g., `m == n` → call built-in `gmres`, `d == 0` → restarted `gmres(m)`)
3. Default parameter assignment (`nargin < k || isempty(param)` idiom)
4. Arnoldi iteration via `modified_gram_schmidt_arnoldi` or `augmented_gram_schmidt_arnoldi`
5. QR factorization via `plane_rotations` → solve least-squares → update `x`
6. Convergence check against relative residual; on failure, compute augmentation vectors and restart

**Test file pattern:** Each `tests/test_*.m` uses MOxUnit's `initTestSuite` / `localfunctions` harness. Tests assert using `assertElementsAlmostEqual`, `assertEqual`, and standard MATLAB `assert`. Test matrices from `data/` are loaded with `load('data/X.mat', 'Problem')`.

## Test data

The `data/` directory holds `.mat` files from the SuiteSparse Matrix Collection. Each file stores a `Problem` struct with fields `A` (sparse matrix) and `b` (right-hand side).

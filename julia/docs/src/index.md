# KrySBAS.jl

**KrySBAS** provides adaptive iterative solvers for sparse linear systems (*Ax = b*) based on Krylov subspaces.

All three solvers share the same output signature:

```julia
x, flag, relresvec, kdvec, time = solver(A, b; kwargs...)
```

| Output | Type | Description |
|--------|------|-------------|
| `x` | `Vector` | Approximate solution |
| `flag` | `Bool` | `true` if converged within `maxit` restarts |
| `relresvec` | `Vector` | Relative residual norm at each restart cycle |
| `kdvec` | `Vector` | Krylov dimension used at each cycle |
| `time` | `Float64` | Elapsed wall-clock time (seconds) |

## Quick start

```julia
using KrySBAS

# LGMRES with m=30 Krylov steps and l=3 error vectors
x, flag, relresvec, kdvec, t = lgmres(A, b; m=30, l=3, tol=1e-8)

# GMRES-E with m=27 and d=3 harmonic Ritz vectors
x, flag, relresvec, kdvec, t = gmres_e(A, b; m=27, d=3, tol=1e-10)

# PD-GMRES with adaptive restart starting at m=20
x, flag, relresvec, kdvec, t = pd_gmres(A, b; m_initial=20, tol=1e-8)
```

## Solvers

| Solver | Description |
|--------|-------------|
| [`gmres_e`](@ref) | Restarted GMRES augmented with harmonic Ritz vectors |
| [`lgmres`](@ref) | Restarted GMRES augmented with error approximation vectors |
| [`pd_gmres`](@ref) | Restarted GMRES with PD-controller restart adaptation |

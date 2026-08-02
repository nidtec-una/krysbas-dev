# KrySBAS.jl

**KrySBAS** provides adaptive iterative solvers for sparse linear systems (*Ax = b*) based on Krylov subspaces.

All solvers share the same output signature:

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

# GMRES-E with m=27 and d=3 harmonic Ritz vectors
x, flag, relresvec, kdvec, t = gmres_e(A, b; m=27, d=3, tol=1e-10)

# GMRES-DR with m=27 Krylov steps and k=3 recycled Ritz vectors
x, flag, relresvec, kdvec, t = gmres_dr(A, b; m=27, k=3, tol=1e-10)

# LGMRES with m=30 Krylov steps and l=3 error vectors
x, flag, relresvec, kdvec, t = lgmres(A, b; m=30, l=3, tol=1e-8)

# PD-GMRES with adaptive restart starting at m=20
x, flag, relresvec, kdvec, t = pd_gmres(A, b; m_initial=20, tol=1e-8)

# SLGMRES-E, switching per cycle between LGMRES-style and GMRES-E-style augmentation
x, flag, relresvec, kdvec, t = slgmres_e(A, b; m=27, l=3, d=3, tol=1e-10)

# A-SLGMRES-E, additionally growing m via a PD law on stagnating cycles
x, flag, relresvec, kdvec, t = a_slgmres_e(A, b; m_initial=27, l=3, d=3, tol=1e-10)
```

## Solvers

Listed chronologically, by the publication year of the method each one implements.

| Solver | Description |
|--------|-------------|
| [`gmres_e`](@ref) | Restarted GMRES augmented with harmonic Ritz vectors ([Morgan, 1995](https://epubs.siam.org/doi/abs/10.1137/S0895479893253975)) |
| [`gmres_dr`](@ref) | Restarted GMRES with deflated (thick) restarting ([Morgan, 2002](https://epubs.siam.org/doi/10.1137/S1064827599364659)) |
| [`lgmres`](@ref) | Restarted GMRES augmented with error approximation vectors ([Baker, Jessup & Manteuffel, 2005](https://epubs.siam.org/doi/abs/10.1137/S0895479803422014)) |
| [`pd_gmres`](@ref) | Restarted GMRES with PD-controller restart adaptation ([Núñez, Schaerer & Bhaya, 2018](https://www.sciencedirect.com/science/article/pii/S037704271830030X)) |
| [`slgmres_e`](@ref) | Restarted GMRES that switches, cycle by cycle, between LGMRES-style and GMRES-E-style augmentation ([Cabral, Schaerer & Bhaya, 2020](https://doi.org/10.1002/nla.2305)) |
| [`a_slgmres_e`](@ref) | `slgmres_e` plus PD-adaptive restart growth on stagnating cycles ([Cabral, Schaerer & Bhaya, 2020](https://doi.org/10.1002/nla.2305)) |

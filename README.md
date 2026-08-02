[![MATLAB tests](https://github.com/nidtec-una/krysbas-dev/actions/workflows/matlab_tests.yaml/badge.svg)](https://github.com/nidtec-una/krysbas-dev/actions/workflows/matlab_tests.yaml)
[![miss_hit](https://github.com/nidtec-una/krysbas-dev/actions/workflows/code_style.yml/badge.svg)](https://github.com/nidtec-una/krysbas-dev/actions/workflows/code_style.yml)
[![Julia tests](https://github.com/nidtec-una/krysbas-dev/actions/workflows/julia_tests.yml/badge.svg)](https://github.com/nidtec-una/krysbas-dev/actions/workflows/julia_tests.yml)
[![Julia style](https://github.com/nidtec-una/krysbas-dev/actions/workflows/julia_style.yml/badge.svg)](https://github.com/nidtec-una/krysbas-dev/actions/workflows/julia_style.yml)
[![codecov](https://codecov.io/gh/nidtec-una/krysbas-dev/graph/badge.svg?token=SRZNZEIBB7)](https://codecov.io/gh/nidtec-una/krysbas-dev)
[![readthedocs](https://img.shields.io/readthedocs/krysbas-dev)](https://krysbas-dev.readthedocs.io/en/latest/?badge=latest)
[![Julia docs](https://github.com/nidtec-una/krysbas-dev/actions/workflows/julia_docs.yml/badge.svg)](https://github.com/nidtec-una/krysbas-dev/actions/workflows/julia_docs.yml)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

# KrySBAS: Krylov Subspace-Based Adaptive Solvers

<p align="center">
  <img src="krysbas_logo.png" alt="KrySBAS logo" width="250"/>
</p>


KrySBAS is a free and open-source toolbox of adaptive iterative solvers for sparse linear systems (*Ax = b*), based on Krylov subspaces. It is available for both **MATLAB/GNU-Octave** and **Julia**.

The toolbox is developed by the [Scientific Computing and Applied Mathematics](https://nidtec.pol.una.py/ccyma/) group at the [NIDTEC](https://nidtec.pol.una.py/) research center of the [Polytechnic Faculty, National University of Asunción, Paraguay](https://www.pol.una.py/).

## Installation

### MATLAB

Clone this repository and add the source directory to your MATLAB path:

```matlab
addpath(genpath('matlab/src'))
```

### Julia

The Julia package lives in the `julia/` subdirectory. Activate and instantiate it once:

```julia
using Pkg
Pkg.activate("julia/")
Pkg.instantiate()
```

Then load the package in your code:

```julia
using KrySBAS
```

## Solvers catalogue

All solvers share the same output signature: `x, flag, relresvec, kdvec, time`. Solvers are listed chronologically, by the publication year of the method each one implements.

### GMRES-E(*m, d*) — [Morgan, 1995](https://epubs.siam.org/doi/abs/10.1137/S0895479893253975)

Restarted GMRES augmented with *d* harmonic Ritz vectors approximating the smallest eigenvalues of the Krylov subspace.

**MATLAB**
```matlab
[x, flag, relresvec, kdvec, time] = gmres_e(A, b, m, d, tol, maxit, xInitial, eigstol)
```

**Julia**
```julia
x, flag, relresvec, kdvec, time = gmres_e(A, b; m=10, d=3, tol=1e-6, maxit=10, x_initial=zeros(n))
```

### GMRES-DR(*m, k*) — [Morgan, 2002](https://epubs.siam.org/doi/10.1137/S1064827599364659)

Restarted GMRES with deflated (thick) restarting: *k* harmonic Ritz vectors are recycled across restart cycles, keeping the subspace dimension fixed at *m* per cycle instead of augmenting it.

**MATLAB**
```matlab
[x, flag, relresvec, kdvec, time, stats] = gmres_dr(A, b, m, k, tol, maxit, xInitial)
```

**Julia**
```julia
x, flag, relresvec, kdvec, time = gmres_dr(A, b; m=10, k=3, tol=1e-6, maxit=10, x_initial=zeros(n))
```

### LGMRES(*m, l*) — [Baker, Jessup & Manteuffel, 2005](https://epubs.siam.org/doi/abs/10.1137/S0895479803422014)

Restarted GMRES augmented with *l* error approximation vectors from prior restart cycles, preserving information from discarded search subspaces.

**MATLAB**
```matlab
[x, flag, relresvec, kdvec, time] = lgmres(A, b, m, l, tol, maxit, xInitial)
```

**Julia**
```julia
x, flag, relresvec, kdvec, time = lgmres(A, b; m=10, l=3, tol=1e-6, maxit=10, x_initial=zeros(n))
```

### PD-GMRES(*m*) — [Núñez, Schaerer & Bhaya, 2018](https://www.sciencedirect.com/science/article/pii/S037704271830030X)

Restarted GMRES with a Proportional-Derivative (PD) controller that automatically adapts the restart parameter *m* each cycle.

**MATLAB**
```matlab
[x, flag, relresvec, kdvec, time] = pd_gmres(A, b, mInitial, mMinMax, mStep, tol, maxit, xInitial, alphaPD)
```

**Julia**
```julia
x, flag, relresvec, kdvec, time = pd_gmres(A, b; m_initial=10, m_min_max=nothing, m_step=1,
                                             tol=1e-6, maxit=10, x_initial=zeros(n), alpha_pd=[-3.0, 5.0])
```

### SLGMRES-E(*m, l, d*) / A-SLGMRES-E(*mⱼ, l, d*) — [Cabral, Schaerer & Bhaya, 2020](https://doi.org/10.1002/nla.2305)

Restarted GMRES that switches, cycle by cycle, between LGMRES-style and GMRES-E-style augmentation based on a convergence-slowdown signal read directly off the residual, at no extra cost. A-SLGMRES-E additionally grows the restart parameter *m* via the same Proportional-Derivative law used by PD-GMRES whenever slowdown is detected.

**MATLAB**
```matlab
[x, flag, relresvec, kdvec, time] = slgmres_e(A, b, m, l, d, epsilonThreshold, tol, maxit, xInitial, eigstol)
[x, flag, relresvec, kdvec, time] = a_slgmres_e(A, b, mInitial, mMinMax, mStep, l, d, epsilonThreshold, alphaPD, tol, maxit, xInitial, eigstol)
```

**Julia**
```julia
x, flag, relresvec, kdvec, time = slgmres_e(A, b; m=10, l=3, d=3, epsilon_threshold=0.01,
                                              tol=1e-6, maxit=10, x_initial=zeros(n))
x, flag, relresvec, kdvec, time = a_slgmres_e(A, b; m_initial=10, m_min_max=nothing, m_step=1, l=3, d=3,
                                                epsilon_threshold=0.01, alpha_pd=[2.0, 0.8],
                                                tol=1e-6, maxit=10, x_initial=zeros(n))
```

## Contributing

If you wish to contribute to KrySBAS, please read the [developer guide](https://github.com/nidtec-una/krysbas-dev/blob/dev_guide/dev_guide.md) before opening a pull request.

## Feature requests and bug reports

For feature requests and bug reports, please create an [issue](https://github.com/nidtec-una/krysbas-dev/issues). For bug reports, please provide a minimal working example that reproduces the error.

[![MATLAB tests](https://github.com/nidtec-una/krysbas-dev/actions/workflows/matlab_tests.yaml/badge.svg)](https://github.com/nidtec-una/krysbas-dev/actions/workflows/matlab_tests.yaml)
[![code_style](https://github.com/nidtec-una/krysbas-dev/actions/workflows/code_style.yml/badge.svg)](https://github.com/nidtec-una/krysbas-dev/actions/workflows/code_style.yml)
[![codecov](https://codecov.io/gh/nidtec-una/krysbas-dev/graph/badge.svg?token=SRZNZEIBB7)](https://codecov.io/gh/nidtec-una/krysbas-dev)
[![docs](https://readthedocs.org/projects/krysbas-dev/badge/?version=latest)](https://krysbas-dev.readthedocs.io/en/latest/?badge=latest)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

# KrySBAS: Krylov Subspace-Based Adaptive Solvers

https://github.com/nidtec-una/krysbas-dev/assets/11464732/c434e9d9-00bf-4446-b035-ed605d86dacd

KrySBAS is a free and open-source MATLAB toolbox containing a collection of adaptive solvers based on Krylov subspaces.  

The toolbox is developed by the [Scientific Computing and Applied Mathematics](https://nidtec.pol.una.py/ccyma/) group at the [NIDTEC](https://nidtec.pol.una.py/) research center of the [Polytechnic Faculty, National University of Asunción, Paraguay](https://www.pol.una.py/).

## Installation

To install KrySBAS, simply clone this repository and add it to your MATLAB path.

## Solvers catalogue

### PD-GMRES(*m*) ([Núñez & Schaerer & Bhaya, 2018](https://www.sciencedirect.com/science/article/pii/S037704271830030X))

Variant of the restarted GMRES that employs a Proportional-Derivative (PD) controller for the automatic selection of the restart parameter *m*.

```Matlab
[x, flag, relresvec, mvec, time] = pd_gmres(A, b, mInitial, mMinMax, mStep, tol, maxit, xInitial, alphaPD)
```

## Contributing

If you wish to contribute to KrySBAS, please read the [developer guide](https://github.com/nidtec-una/krysbas-dev/blob/dev_guide/dev_guide.md) before opening a pull request.

## Feature requests and bug reports

For future requests and bug reports, please create an [issue](https://github.com/nidtec-una/krysbas-dev/issues). In the latter case, we kindly ask you to provide a MWE that reproduces the error.

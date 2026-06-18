# Minimal example: solve Ax = b with GMRES using Krylov.jl
#
# Run from VS Code by opening this file and pressing Alt+Enter (line/block)
# or Ctrl+F5 to run the whole file.
#
# No prior setup needed — a temporary Julia environment is created
# automatically and Krylov.jl is installed into it on first run.

import Pkg
Pkg.activate(; temp=true)
Pkg.add("Krylov"; io=devnull)

using Krylov
using LinearAlgebra

# ── Build a 1-D Poisson problem: -u'' = f on [0,1] with Dirichlet BCs ────────
# Discretised with n interior points → tridiagonal system
n = 100
h = 1.0 / (n + 1)

diag_main = fill(2.0 / h^2, n)
diag_off  = fill(-1.0 / h^2, n - 1)
A = diagm(0 => diag_main, 1 => diag_off, -1 => diag_off)

# Right-hand side: f(x) = π² sin(πx), exact solution u(x) = sin(πx)
x_grid = [(i * h) for i in 1:n]
b = π^2 .* sin.(π .* x_grid)

# ── Solve with GMRES ──────────────────────────────────────────────────────────
x, stats = Krylov.gmres(A, b; memory=20, rtol=1e-8, verbose=0)

# ── Report ────────────────────────────────────────────────────────────────────
u_exact = sin.(π .* x_grid)
err     = norm(x - u_exact, Inf)

println("Converged : ", stats.solved)
println("Iterations: ", stats.niter)
println("Residual  : ", norm(b - A * x) / norm(b))
println("Max error vs exact solution: ", round(err; sigdigits=4))

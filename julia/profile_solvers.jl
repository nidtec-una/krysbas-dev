#!/usr/bin/env julia
# Quick allocation and timing profile for the two solvers on sherman5.
# Run from julia/ with: julia --project=. profile_solvers.jl

using KrySBAS, MAT, LinearAlgebra, SparseArrays, Printf

DATA = joinpath(@__DIR__, "..", "data")
f = matopen(joinpath(DATA, "sherman5.mat")); P = read(f,"Problem"); close(f)
A = P["A"]; b = vec(P["b"])

TOL = 1e-9; MAXIT = 1000

println("=== lgmres (m=100, l=3) on sherman5 ===")
lgmres(A, b; m=100, l=3, tol=TOL, maxit=MAXIT)          # warmup
alloc = @allocated lgmres(A, b; m=100, l=3, tol=TOL, maxit=MAXIT)
t     = @elapsed    lgmres(A, b; m=100, l=3, tol=TOL, maxit=MAXIT)
@printf "  time: %.3f s   allocations: %.1f MB\n" t alloc/1e6

println()
println("=== pd_gmres (m_initial=30, m_step=3) on sherman5 ===")
pd_gmres(A, b; m_initial=30, m_step=3, tol=TOL, maxit=MAXIT)  # warmup
alloc = @allocated pd_gmres(A, b; m_initial=30, m_step=3, tol=TOL, maxit=MAXIT)
t     = @elapsed    pd_gmres(A, b; m_initial=30, m_step=3, tol=TOL, maxit=MAXIT)
@printf "  time: %.3f s   allocations: %.1f MB\n" t alloc/1e6

println()
println("=== Inner kernels on a synthetic n=3312, s=30 problem ===")
import KrySBAS: augmented_gram_schmidt_arnoldi, modified_gram_schmidt_arnoldi, plane_rotations

n = 3312; m = 27; k = 3
v1 = randn(n); v1 ./= norm(v1)
appendV = randn(n, k)
for col in eachcol(appendV); col ./= norm(col); end

# warmup
augmented_gram_schmidt_arnoldi(A, v1, m, appendV)
modified_gram_schmidt_arnoldi(A, v1, m)

a_alloc = @allocated augmented_gram_schmidt_arnoldi(A, v1, m, appendV)
a_time  = @elapsed   augmented_gram_schmidt_arnoldi(A, v1, m, appendV)
@printf "  augmented_gram_schmidt_arnoldi:  %.4f s  %.1f MB\n" a_time a_alloc/1e6

m_alloc = @allocated modified_gram_schmidt_arnoldi(A, v1, m)
m_time  = @elapsed   modified_gram_schmidt_arnoldi(A, v1, m)
@printf "  modified_gram_schmidt_arnoldi:   %.4f s  %.1f MB\n" m_time m_alloc/1e6

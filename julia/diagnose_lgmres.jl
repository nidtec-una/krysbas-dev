#!/usr/bin/env julia
# Diagnostic: compare lgmres residual history Julia vs Octave on all three matrices.
# Run from julia/ with: julia --project=. diagnose_lgmres.jl

using KrySBAS, MAT, LinearAlgebra, SparseArrays, Printf

const TOL   = 1e-9
const MAXIT = 1000
const M     = 27
const L     = 3
const DATA  = joinpath(@__DIR__, "..", "data")
const OCT   = joinpath(@__DIR__, "..", "matlab")

# ---- Run Octave and collect per-matrix residual vectors ----
function run_octave_diag(mat_name)
    script = tempname() * ".m"
    octave_bin = something(Sys.which("octave"), "octave")
    open(script, "w") do io
        println(io, "warning('off','all');")
        println(io, "addpath(genpath('$(OCT)/src'));")
        println(io, "load('$(DATA)/$(mat_name).mat','Problem');")
        println(io, "A=Problem.A; b=Problem.b;")
        println(io, "[~,~,rrv,~,~]=lgmres(A,b,$(M),$(L),$(TOL),$(MAXIT));")
        println(io, "for k=1:length(rrv); fprintf('RRV %.16e\\n',rrv(k)); end")
    end
    cmd = pipeline(`$octave_bin --norc --no-gui $script`; stderr=devnull)
    raw = try read(cmd, String) catch; return nothing end
    rm(script; force=true)
    vals = Float64[]
    for line in split(raw, '\n')
        startswith(line, "RRV ") || continue
        push!(vals, parse(Float64, strip(line[5:end])))
    end
    return vals
end

# ---- Print comparison table ----
function compare(mat_name)
    println("\n" * "="^70)
    println("  lgmres on $mat_name  (m=$M, l=$L, tol=$TOL, maxit=$MAXIT)")
    println("="^70)

    f   = matopen(joinpath(DATA, "$mat_name.mat"))
    P   = read(f, "Problem"); close(f)
    A   = P["A"]; b = vec(P["b"])

    _, _, jl_rrv, _, _ = lgmres(A, b; m=M, l=L, tol=TOL, maxit=MAXIT)
    oct_rrv = run_octave_diag(mat_name)

    if oct_rrv === nothing
        println("  Octave run failed.")
        return
    end

    n_jl  = length(jl_rrv)
    n_oct = length(oct_rrv)
    @printf "  Julia:  %d entries, final relres = %.4e\n" n_jl  jl_rrv[end]
    @printf "  Octave: %d entries, final relres = %.4e\n" n_oct oct_rrv[end]

    n_common = min(n_jl, n_oct)
    diffs    = [abs(jl_rrv[k] - oct_rrv[k]) for k in 1:n_common]
    max_diff = maximum(diffs)
    max_idx  = argmax(diffs)
    @printf "  Max |Julia - Octave| = %.4e at cycle %d\n" max_diff max_idx

    println()
    println("  cycle     Julia relres    Octave relres    |diff|")
    println("  " * "-"^58)
    show_idx = unique(vcat(1:min(5,n_common), max(1,n_common-4):n_common, max_idx))
    sort!(show_idx)
    prev = 0
    for k in show_idx
        if k > prev + 1
            println("  ...")
        end
        @printf "  %5d   %14.6e   %14.6e   %10.2e\n" k jl_rrv[k] oct_rrv[k] diffs[k]
        prev = k
    end

    # For non-converging cases: find first stagnation cycle in each
    if jl_rrv[end] >= TOL
        jl_stag = findfirst(k -> abs(jl_rrv[k] - jl_rrv[k-1]) < 1e-12,
                            2:n_jl)
        oct_stag = findfirst(k -> abs(oct_rrv[k] - oct_rrv[k-1]) < 1e-12,
                             2:n_oct)
        jl_stag  = isnothing(jl_stag)  ? "none" : string(jl_stag  + 1)
        oct_stag = isnothing(oct_stag) ? "none" : string(oct_stag + 1)
        println()
        println("  First stagnation cycle (Δrelres < 1e-12):")
        println("    Julia: $jl_stag    Octave: $oct_stag")
    end
end

for mat in ["sherman1", "sherman4", "sherman5"]
    isfile(joinpath(DATA, "$mat.mat")) && compare(mat)
end
println()

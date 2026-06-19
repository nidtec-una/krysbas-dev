#!/usr/bin/env julia
# KrySBAS performance comparison: Julia vs Octave
#
# Usage (from the repo root):
#   julia --project=julia julia/benchmark.jl
#
# Usage (from julia/):
#   julia --project=. benchmark.jl

using KrySBAS, MAT, LinearAlgebra, SparseArrays, Printf

const MATRICES = ["sherman1", "sherman4", "sherman5"]
const TOL = 1e-9
const MAXIT = 1000
const N_RUNS = 3   # timed runs (after 1 warmup) — minimum is reported
const DATA_DIR = joinpath(@__DIR__, "..", "data")

# Solver dispatch table with fixed parameters (gmres_e excluded: SM vs LM Ritz mismatch)
const SOLVERS = [
    ("lgmres", (A, b) -> lgmres(A, b; m = 100, l = 3, tol = TOL, maxit = MAXIT)),
    (
        "pd_gmres",
        (A, b) -> pd_gmres(A, b; m_initial = 30, m_step = 3, tol = TOL, maxit = MAXIT),
    ),
]

function load_mat(name)
    f = matopen(joinpath(DATA_DIR, "$name.mat"))
    P = read(f, "Problem");
    close(f)
    return P["A"], vec(P["b"])
end

# Returns (min_time_s, iters, relres, flag) or nothing on failure
function bench_solver(fn, A, b)
    try
        fn(A, b)
    catch
        ;
        return nothing
    end   # warmup: triggers JIT compilation
    best = Inf
    local iters, relres, flag
    for _ = 1:N_RUNS
        t0 = time_ns()
        _, fl, rrv, _, _ = fn(A, b)
        elapsed = (time_ns() - t0) / 1e9
        best = min(best, elapsed)
        iters = length(rrv) - 1
        relres = rrv[end]
        flag = fl
    end
    return best, iters, relres, flag
end

# ---- Run Octave benchmark ----
function run_octave()
    octave_script = joinpath(@__DIR__, "..", "matlab", "benchmark.m")
    if !isfile(octave_script)
        @warn "matlab/benchmark.m not found — skipping Octave."
        return nothing
    end
    octave_bin = something(Sys.which("octave"), "octave")
    cmd = pipeline(`$octave_bin --norc --no-gui $octave_script`; stderr = devnull)
    output = try
        read(cmd, String)
    catch e
        @warn "Octave run failed: $e"
        return nothing
    end
    # Parse RESULT lines: RESULT,solver,matrix,n,nnz,time,iters,relres
    results = Dict()
    for line in split(output, '\n')
        startswith(line, "RESULT,") || continue
        parts = split(strip(line), ',')
        length(parts) == 8 || continue
        _, solver, mat, n, nnz_a, t, iters, relres = parts
        key = (solver, mat)
        results[key] = (
            t = parse(Float64, t),
            iters = parse(Int, iters),
            relres = parse(Float64, relres),
            n = parse(Int, n),
            nnz = parse(Int, nnz_a),
        )
    end
    return results
end

# ---- Formatting helpers ----
hr(w) = "─" ^ w
hdr(w) = "═" ^ w

function print_table(mat_name, n, nnz_a, julia_res, oct_res)
    println()
    @printf "  Matrix: %-10s  n = %-6d  nnz = %d\n" mat_name n nnz_a
    @printf "  %-12s %12s %12s %9s %9s %9s %9s\n" "Solver" "Julia t(s)" "Octave t(s)" "Speedup" "J iters" "O iters" "Rel.res."
    println("  " * hr(85))
    for (sname, _) in SOLVERS
        has_jl = haskey(julia_res, sname)
        has_oct = oct_res !== nothing && haskey(oct_res, (sname, mat_name))
        jt = has_jl ? julia_res[sname][1] : NaN
        ji = has_jl ? julia_res[sname][2] : -1
        jr = has_jl ? julia_res[sname][3] : NaN
        jt_str = has_jl ? @sprintf("%12.4f", jt) : @sprintf("%12s", "ERR")
        ji_str = has_jl ? @sprintf("%9d", ji) : @sprintf("%9s", "ERR")
        jr_str = has_jl ? @sprintf("%9.2e", jr) : @sprintf("%9s", "ERR")
        if has_oct
            ot = oct_res[(sname, mat_name)].t
            oi = oct_res[(sname, mat_name)].iters
            ot_str = @sprintf "%12.4f" ot
            oi_str = @sprintf "%9d" oi
            sp_str = has_jl ? @sprintf("%8.1fx", ot / jt) : @sprintf("%9s", "n/a")
        else
            ot_str = @sprintf "%12s" "n/a"
            oi_str = @sprintf "%9s" "n/a"
            sp_str = @sprintf "%9s" "n/a"
        end
        @printf "  %-12s %s %s %s %s %s %s\n" sname jt_str ot_str sp_str ji_str oi_str jr_str
    end
end

# ---- Main ----
println()
println(hdr(70))
println("  KrySBAS Performance Comparison  ·  Julia vs Octave")
println("  Parameters: lgmres m=100 l=3  |  pd_gmres m_initial=30 m_step=3  |  tol=1e-9")
println("  Julia timing: warmup + $(N_RUNS) runs, minimum reported")
println("  Octave timing: $(N_RUNS) runs, minimum reported")
println(hdr(70))

print("\nRunning Octave benchmark...")
flush(stdout)
oct_results = run_octave()
println(oct_results === nothing ? " (skipped)" : " done.")

# Pre-load matrices once (for warmup across all solvers)
for mat_name in MATRICES
    isfile(joinpath(DATA_DIR, "$mat_name.mat")) || continue
    A, b = load_mat(mat_name)
    n, nnz_a = size(A, 1), nnz(A)

    julia_res = Dict()
    for (sname, fn) in SOLVERS
        print("Benchmarking Julia $sname on $mat_name...")
        flush(stdout)
        result = bench_solver(fn, A, b)
        if result === nothing
            println(" ERROR (NaN/Inf — skipped)")
        else
            t, iters, relres, flag = result
            julia_res[sname] = (t, iters, relres, flag)
            println(
                @sprintf " %.4fs, %d iters, relres=%.2e%s" t iters relres (
                    flag ? "" : " [did not converge]"
                )
            )
        end
    end

    print_table(mat_name, n, nnz_a, julia_res, oct_results)
end

println()
println(hdr(70))

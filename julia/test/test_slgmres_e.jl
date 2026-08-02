@testset "slgmres_e" begin

    @testset "input validation" begin
        @test_throws ArgumentError slgmres_e(zeros(0, 0), Float64[])
        @test_throws ArgumentError slgmres_e([1.0 2.0; 3.0 4.0; 5.0 6.0], [1.0; 1.0; 1.0])
        @test_throws ArgumentError slgmres_e(Matrix{Float64}(I, 3, 3), Float64[])
        @test_throws ArgumentError slgmres_e(Matrix{Float64}(I, 3, 3), [1.0; 1.0])
        @test_throws ArgumentError slgmres_e(Matrix{Float64}(I, 3, 3), ones(3); m = 4)
        @test_throws ArgumentError slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m = 2,
            l = 0,
        )
        @test_throws ArgumentError slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m = 2,
            l = 1,
            d = 0,
        )
        @test_throws ArgumentError slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m = 2,
            l = 1,
            d = 1,
            epsilon_threshold = 1.5,
        )
        @test_throws ArgumentError slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            x_initial = ones(2),
        )
    end

    @testset "m == n dispatches to full unrestarted GMRES" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = slgmres_e(A, b; m = 3)
        @test x ≈ ones(3) atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 3)
        @test t > 0
    end

    @testset "default parameters: identity matrix" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, _, _, _ = slgmres_e(A, b)
        @test x ≈ ones(3) atol=1e-10
        @test flag
    end

    @testset "identity n=3 m=2 l=1 d=1: happy breakdown on cycle 1" begin
        # A*v1 is parallel to v1 for A = I, so the Arnoldi loop hits a
        # happy breakdown after a single step regardless of m.
        A = Matrix{Float64}(I, 3, 3)
        b = [2.0; 3.0; 4.0]
        x, flag, relresvec, kdvec, t =
            slgmres_e(A, b; m = 2, l = 1, d = 1, tol = 1e-9, maxit = 100)
        @test x ≈ [2.0; 3.0; 4.0] atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test kdvec[1] == 1
        @test t > 0
    end

    @testset "Embree 3x3 toy example" begin
        # slgmres_e.m (Octave) converges on this case in 2 cycles, but
        # that is numerical luck rather than a target to match here. At
        # n=3, m=2, l=1, d=1, both augmentation branches (LGMRES-style
        # z_mat and GMRES-E-style dy) propose a single extra direction
        # that is, up to machine precision, already in the span of the
        # m=2 Krylov basis (h/w_norm ~ 2e-16 in
        # augmented_gram_schmidt_arnoldi's near-breakdown check -- a
        # genuine, not borderline, redundancy). augmented_gram_schmidt
        # _arnoldi.jl correctly detects this and declines to augment,
        # which permanently caps the subspace at s=m=2 and stalls the
        # solver (same structural-stall class as gmres_dr's own embree3
        # test). The unguarded MATLAB/Octave Arnoldi has no such check,
        # so it proceeds with a numerically redundant direction anyway
        # and, for this specific exactly-representable 3x3 system,
        # happens to still land on the exact answer. This test documents
        # the known, expected stall in Julia rather than asserting
        # convergence.
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "embree3.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        x, flag, relresvec, kdvec, t =
            slgmres_e(A, b; m = 2, l = 1, d = 1, tol = 1e-6, maxit = 100)
        @test !flag
        @test all(==(2), kdvec)
        @test relresvec[end] ≈ 0.3764959840191753 atol=1e-8
        @test t > 0
    end

    @testset "Sherman1 matches slgmres_e.m exactly" begin
        # Regression test against the validated MATLAB slgmres_e.m: both
        # converge in exactly 38 cycles to matching precision on this
        # problem, confirming the port's core switching/augmentation logic
        # (not just the standalone LGMRES/GMRES-E pieces it is built from).
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman1.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        x, flag, relresvec, kdvec, t =
            slgmres_e(A, b; m = 27, l = 3, d = 3, tol = 1e-12, maxit = 1000)
        @test flag
        @test length(relresvec) - 1 == 38
        @test relresvec[end] < 1e-12
        @test t > 0
    end

    @testset "Sherman4" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman4.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, _, t =
            slgmres_e(A, b; m = 27, l = 3, d = 3, tol = 1e-12, maxit = 1000)
        @test flag
        @test relresvec[end] < 1e-12
        @test t > 0
    end

    @testset "Sherman5 converges (tight m=28,l=2,d=2 budget)" begin
        # m - k_min = 1 here (only one fresh direction once l/d saturate at
        # their minimum), the same flavor of tight-budget case that showed
        # cross-language numerical sensitivity in the harmonic-Ritz
        # eigenproblem for gmres_dr on embree3 earlier in this port. Both
        # this port and slgmres_e.m converge correctly here, just on a
        # different cycle schedule (not asserted exactly, unlike Sherman1).
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman5.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, _, t = slgmres_e(
            A,
            b;
            m = 28,
            l = 2,
            d = 2,
            epsilon_threshold = 0.01,
            tol = 1e-9,
            maxit = 1000,
        )
        @test flag
        @test relresvec[end] < 1e-9
        @test t > 0
    end

end

@testset "gmres_dr" begin

    @testset "input validation" begin
        @test_throws ArgumentError gmres_dr(zeros(0, 0), Float64[])
        @test_throws ArgumentError gmres_dr([1.0 2.0; 3.0 4.0; 5.0 6.0], [1.0; 1.0; 1.0])
        @test_throws ArgumentError gmres_dr(Matrix{Float64}(I, 3, 3), Float64[])
        @test_throws ArgumentError gmres_dr(Matrix{Float64}(I, 3, 3), [1.0; 1.0])
        @test_throws ArgumentError gmres_dr(Matrix{Float64}(I, 3, 3), ones(3); m = 4)
        @test_throws ArgumentError gmres_dr(Matrix{Float64}(I, 3, 3), ones(3); m = 2, k = 2)
        @test_throws ArgumentError gmres_dr(Matrix{Float64}(I, 3, 3), ones(3); m = 2, k = 3)
        @test_throws ArgumentError gmres_dr(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            x_initial = ones(2),
        )
    end

    @testset "m == n dispatches to full unrestarted GMRES" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = gmres_dr(A, b; m = 3)
        @test x ≈ ones(3) atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 3)
        @test t > 0
    end

    @testset "k == 0 dispatches to restarted GMRES(m)" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = gmres_dr(A, b; m = 2, k = 0)
        @test x ≈ ones(3) atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 2)
        @test t > 0
    end

    @testset "default parameters: identity matrix" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, _, _, _ = gmres_dr(A, b)
        @test x ≈ ones(3) atol=1e-10
        @test flag
    end

    @testset "identity n=3 m=2 k=1: happy breakdown on iteration 1" begin
        A = Matrix{Float64}(I, 3, 3)
        b = [2.0; 3.0; 4.0]
        x, flag, relresvec, kdvec, t = gmres_dr(A, b; m = 2, k = 1, tol = 1e-9, maxit = 100)
        @test x ≈ [2.0; 3.0; 4.0] atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test kdvec[1] == 1
        @test t > 0
    end

    @testset "diagonal matrix n=10 m=4 k=2" begin
        n = 10
        A = diagm(1.0:n)
        b = ones(n)
        x, flag, _, _, t = gmres_dr(A, b; m = 4, k = 2, tol = 1e-10, maxit = 200)
        @test x ≈ (1.0 ./ (1:n)) atol=1e-8
        @test flag
        @test t > 0
    end

    @testset "Embree 3x3 toy example: documented stall (m-k=1 too tight)" begin
        # GMRES-DR(m=2, k=1) is the only non-trivial (m, k) combination for a
        # 3x3 system. This matrix's harmonic Ritz pencil at that size is a
        # genuine complex-conjugate pair, so recycling it correctly forces
        # keep = k+1 = m, leaving no budget for a fresh Arnoldi direction in
        # any later cycle: GMRES-DR stalls on this specific problem by
        # mathematical necessity (see matlab/tests/test_gmres_dr.m for the
        # matching MATLAB documentation of the same behaviour).
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "embree3.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, _, t = gmres_dr(A, b; m = 2, k = 1, tol = 1e-6, maxit = 100)
        @test !flag
        @test relresvec[end] ≈ 0.4629 rtol = 1e-3
        @test t > 0
    end

    @testset "Sherman1" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman1.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, kdvec, t =
            gmres_dr(A, b; m = 27, k = 3, tol = 1e-12, maxit = 1000)
        @test flag
        @test relresvec[end] < 1e-12
        @test all(kdvec .<= 27)
        @test t > 0
    end

    @testset "Sherman4" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman4.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, _, t = gmres_dr(A, b; m = 27, k = 3, tol = 1e-12, maxit = 1000)
        @test flag
        @test relresvec[end] < 1e-12
        @test t > 0
    end

end

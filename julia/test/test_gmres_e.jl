@testset "gmres_e" begin

    @testset "input validation" begin
        @test_throws ArgumentError gmres_e(zeros(0, 0), Float64[])
        @test_throws ArgumentError gmres_e([1.0 2.0; 3.0 4.0; 5.0 6.0], [1.0; 1.0; 1.0])
        @test_throws ArgumentError gmres_e(Matrix{Float64}(I, 3, 3), Float64[])
        @test_throws ArgumentError gmres_e(Matrix{Float64}(I, 3, 3), [1.0; 1.0])
        @test_throws ArgumentError gmres_e(Matrix{Float64}(I, 3, 3), ones(3); m = 4)
        @test_throws ArgumentError gmres_e(Matrix{Float64}(I, 3, 3), ones(3); m = 2, d = 3)
        @test_throws ArgumentError gmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            x_initial = ones(2),
        )
    end

    @testset "default m uses min(n, 10) → dispatches to full GMRES for small A" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = gmres_e(A, b)
        @test x ≈ ones(3) atol=1e-10
        @test flag
        @test t > 0
    end

    @testset "m == n dispatches to full unrestarted GMRES" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = gmres_e(A, b; m = 3, d = 0)
        @test x ≈ ones(3) atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 3)
        @test t > 0
    end

    @testset "d == 0 dispatches to restarted GMRES(m)" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = gmres_e(A, b; m = 2, d = 0)
        @test x ≈ ones(3) atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 2)
        @test t > 0
    end

    @testset "identity m=1 d=1: happy breakdown on iteration 1" begin
        A = Matrix{Float64}(I, 3, 3)
        b = [2.0; 3.0; 4.0]
        x, flag, relresvec, kdvec, t = gmres_e(A, b; m = 1, d = 1, tol = 1e-9, maxit = 100)
        @test x ≈ [2.0; 3.0; 4.0] atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 1)
        @test t > 0
    end

    @testset "issue 68: identity m=2 d=1 happy breakdown" begin
        A = Matrix{Float64}(I, 3, 3)
        b = [2.0; 3.0; 4.0]
        x, flag, relresvec, kdvec, t = gmres_e(A, b; m = 2, d = 1, tol = 1e-9, maxit = 100)
        @test x ≈ [2.0; 3.0; 4.0] atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 1)
        @test t > 0
    end

    @testset "Embree 3x3 toy example" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "embree3.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        x, flag, _, _, t = gmres_e(A, b; m = 2, d = 1, tol = 1e-6, maxit = 100)
        @test x ≈ [8.0; -7.0; 1.0] atol=1e-4
        @test flag
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
            gmres_e(A, b; m = 27, d = 3, tol = 1e-12, maxit = 1000)
        @test flag
        @test relresvec[end] < 1e-12
        @test all(kdvec .== 30)
        @test t > 0
    end

    @testset "Sherman4" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman4.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, kdvec, t =
            gmres_e(A, b; m = 27, d = 3, tol = 1e-12, maxit = 1000)
        @test flag
        @test relresvec[end] < 1e-12
        @test all(kdvec .== 30)
        @test t > 0
    end

    @testset "Sherman5" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman5.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, kdvec, t =
            gmres_e(A, b; m = 27, d = 3, tol = 1e-12, maxit = 1000)
        @test flag
        @test relresvec[end] < 1e-12
        @test t > 0
    end

end

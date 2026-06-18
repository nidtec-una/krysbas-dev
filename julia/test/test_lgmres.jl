@testset "lgmres" begin

    @testset "input validation" begin
        @test_throws ArgumentError lgmres(zeros(0, 0), Float64[])
        @test_throws ArgumentError lgmres([1.0 2.0; 3.0 4.0; 5.0 6.0], [1.0; 1.0; 1.0])
        @test_throws ArgumentError lgmres(Matrix{Float64}(I, 3, 3), Float64[])
        @test_throws ArgumentError lgmres(Matrix{Float64}(I, 3, 3), [1.0; 1.0])
        @test_throws ArgumentError lgmres(Matrix{Float64}(I, 3, 3), ones(3); m=4)
        @test_throws ArgumentError lgmres(Matrix{Float64}(I, 3, 3), ones(3);
            x_initial=ones(2))
    end

    @testset "default m dispatches to full GMRES for small A" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = lgmres(A, b)
        @test x ≈ ones(3) atol=1e-10
        @test flag
        @test t > 0
    end

    @testset "m == n dispatches to full unrestarted GMRES" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = lgmres(A, b; m=3, l=0)
        @test x ≈ ones(3) atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 3)
        @test t > 0
    end

    @testset "l == 0 dispatches to restarted GMRES(m)" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = lgmres(A, b; m=2, l=0)
        @test x ≈ ones(3) atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 2)
        @test t > 0
    end

    @testset "identity n=100 m=27 l=3: converges in one cycle" begin
        n = 100
        A = Matrix{Float64}(I, n, n)
        b = ones(n)
        x, flag, relresvec, kdvec, t = lgmres(A, b; m=27, l=3, tol=1e-6, maxit=100)
        @test x ≈ ones(n) atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 30)
        @test t > 0
    end

    @testset "identity n=3 m=2 l=1: happy breakdown on iteration 1" begin
        A = Matrix{Float64}(I, 3, 3)
        b = [2.0; 3.0; 4.0]
        x, flag, relresvec, kdvec, t = lgmres(A, b; m=2, l=1, tol=1e-6, maxit=100)
        @test x ≈ [2.0; 3.0; 4.0] atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(kdvec .== 3)
        @test t > 0
    end

    @testset "Embree 3x3 toy example" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "embree3.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        x, flag, _, _, t = lgmres(A, b; m=2, l=1, tol=1e-6, maxit=100)
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

        _, flag, relresvec, _, t = lgmres(A, b; m=27, l=3, tol=1e-12, maxit=1000)
        @test flag
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

        _, flag, relresvec, _, t = lgmres(A, b; m=27, l=3, tol=1e-12, maxit=1000)
        @test flag
        @test relresvec[end] < 1e-12
        @test t > 0
    end

end

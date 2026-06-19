@testset "pd_gmres" begin

    @testset "input validation" begin
        @test_throws ArgumentError pd_gmres(zeros(0, 0), Float64[])
        @test_throws ArgumentError pd_gmres([1.0 2.0; 3.0 4.0; 5.0 6.0], [1.0; 1.0; 1.0])
        @test_throws ArgumentError pd_gmres(Matrix{Float64}(I, 3, 3), Float64[])
        @test_throws ArgumentError pd_gmres(Matrix{Float64}(I, 3, 3), [1.0; 1.0])
        @test_throws ArgumentError pd_gmres(
            Matrix{Float64}(I, 4, 4),
            ones(4);
            m_initial = -1,
        )
        @test_throws ArgumentError pd_gmres(
            Matrix{Float64}(I, 4, 4),
            ones(4);
            m_initial = 5,
        )
        @test_throws ArgumentError pd_gmres(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 1,
            m_min_max = [0, 2],
        )
        @test_throws ArgumentError pd_gmres(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 1,
            m_min_max = [1, 4],
        )
        @test_throws ArgumentError pd_gmres(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 1,
            m_min_max = [1, 1],
        )
        @test_throws ArgumentError pd_gmres(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 1,
            m_min_max = [2, 3],
        )
        @test_throws ArgumentError pd_gmres(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 1,
            m_step = 0,
        )
        @test_throws ArgumentError pd_gmres(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            x_initial = ones(2),
        )
    end

    @testset "m_initial=0 (unset) or m_initial=n → unrestarted, kdvec=[NaN]" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        # m_initial not given
        _, _, _, kdvec, _ = pd_gmres(A, b)
        @test all(isnan, kdvec)
        # m_initial == n
        _, _, _, kdvec2, _ = pd_gmres(A, b; m_initial = 3)
        @test all(isnan, kdvec2)
    end

    @testset "m_initial provided → restarted, kdvec is Int vector" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        _, _, _, kdvec, _ = pd_gmres(A, b; m_initial = 1)
        @test !all(isnan, kdvec)
    end

    @testset "unrestarted identity matrix" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = pd_gmres(A, b)
        @test x ≈ ones(3) atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test all(isnan, kdvec)
        @test t > 0
    end

    @testset "restarted identity m_initial=1: happy breakdown, converges" begin
        A = Matrix{Float64}(I, 3, 3)
        b = [2.0; 3.0; 4.0]
        x, flag, relresvec, kdvec, t = pd_gmres(A, b; m_initial = 1, tol = 1e-9, maxit = 10)
        @test x ≈ [2.0; 3.0; 4.0] atol=1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol=1e-14
        @test kdvec == [1; 1]
        @test t > 0
    end

    @testset "Embree 3x3, m_initial=1" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "embree3.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        x, flag, relresvec, kdvec, _ = pd_gmres(A, b; m_initial = 1, tol = 1e-9, maxit = 20)
        @test x ≈ [8.0; -7.0; 1.0] atol=1e-4
        @test flag
        @test kdvec == [1; 1; 1; 2]
        @test relresvec ≈ [
            1.0;
            0.925820099772551;
            0.654653670707977;
            0.0
        ] atol=1e-8
    end

    @testset "Embree 3x3, m_initial=2: PD rule escapes stagnation" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "embree3.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        x, flag, relresvec, kdvec, _ = pd_gmres(A, b; m_initial = 2, tol = 1e-9, maxit = 20)
        @test x ≈ [8.0; -7.0; 1.0] atol=1e-4
        @test flag
        @test kdvec == [2; 2; 2; 3]
        @test relresvec[1] ≈ 1.0
        @test relresvec[end] < 1e-9
    end

    @testset "Sherman1" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman1.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, _, _ =
            pd_gmres(A, b; m_initial = 30, m_step = 3, tol = 1e-9, maxit = 1000)
        @test flag
        @test relresvec[end] < 1e-9
    end

    @testset "Sherman4" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman4.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, _, _ =
            pd_gmres(A, b; m_initial = 30, m_step = 3, tol = 1e-9, maxit = 1000)
        @test flag
        @test relresvec[end] < 1e-9
    end

end

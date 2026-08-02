@testset "a_slgmres_e" begin

    @testset "input validation" begin
        @test_throws ArgumentError a_slgmres_e(zeros(0, 0), Float64[])
        @test_throws ArgumentError a_slgmres_e([1.0 2.0; 3.0 4.0; 5.0 6.0], [1.0; 1.0; 1.0])
        @test_throws ArgumentError a_slgmres_e(Matrix{Float64}(I, 3, 3), Float64[])
        @test_throws ArgumentError a_slgmres_e(Matrix{Float64}(I, 3, 3), [1.0; 1.0])
        @test_throws ArgumentError a_slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 4,
        )
        @test_throws ArgumentError a_slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 2,
            m_min_max = [2, 1],
        )
        @test_throws ArgumentError a_slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 2,
            m_step = 0,
        )
        @test_throws ArgumentError a_slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 2,
            l = 0,
        )
        @test_throws ArgumentError a_slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 2,
            l = 1,
            d = 0,
        )
        @test_throws ArgumentError a_slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            m_initial = 2,
            l = 1,
            d = 1,
            epsilon_threshold = 1.5,
        )
        @test_throws ArgumentError a_slgmres_e(
            Matrix{Float64}(I, 3, 3),
            ones(3);
            x_initial = ones(2),
        )
    end

    @testset "m_initial == n dispatches to full unrestarted GMRES" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, relresvec, kdvec, t = a_slgmres_e(A, b; m_initial = 3)
        @test x ≈ ones(3) atol = 1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol = 1e-14
        @test all(kdvec .== 3)
        @test t > 0
    end

    @testset "default parameters: identity matrix" begin
        A = Matrix{Float64}(I, 3, 3)
        b = ones(3)
        x, flag, _, _, _ = a_slgmres_e(A, b)
        @test x ≈ ones(3) atol = 1e-10
        @test flag
    end

    @testset "identity n=3 m_initial=2 l=1 d=1: happy breakdown on cycle 1" begin
        # A*v1 is parallel to v1 for A = I, so the Arnoldi loop hits a
        # happy breakdown after a single step regardless of m.
        A = Matrix{Float64}(I, 3, 3)
        b = [2.0; 3.0; 4.0]
        x, flag, relresvec, kdvec, t =
            a_slgmres_e(A, b; m_initial = 2, l = 1, d = 1, tol = 1e-9, maxit = 100)
        @test x ≈ [2.0; 3.0; 4.0] atol = 1e-10
        @test flag
        @test relresvec ≈ [1.0; 0.0] atol = 1e-14
        @test kdvec[1] == 1
        @test t > 0
    end

    @testset "Embree 3x3 toy example" begin
        # Same structural stall as slgmres_e.jl's own embree3 test (see
        # that file's comment for the full explanation): at n=3,
        # m_initial=2, l=1, d=1, both augmentation branches propose an
        # extra direction that is, up to machine precision, already in
        # the span of the m=2 Krylov basis, so
        # augmented_gram_schmidt_arnoldi's near-breakdown check correctly
        # declines to augment. Here that also means the stagnating branch
        # never gets a chance to grow m past m_initial (pd_rule is only
        # invoked on stagnating cycles, and the solver never leaves the
        # non-stagnating branch once relresvec stops decreasing quickly
        # enough to keep re-triggering it -- see the frozen kdvec below).
        # MATLAB/Octave's unguarded Arnoldi has no such check and,
        # numerically lucky on this exact tiny system, converges in 2
        # cycles.
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "embree3.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        x, flag, relresvec, kdvec, t =
            a_slgmres_e(A, b; m_initial = 2, l = 1, d = 1, tol = 1e-6, maxit = 100)
        @test !flag
        @test all(==(2), kdvec)
        @test t > 0
    end

    @testset "Sherman1" begin
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman1.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, _, t =
            a_slgmres_e(A, b; m_initial = 27, l = 3, d = 3, tol = 1e-12, maxit = 1000)
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

        _, flag, relresvec, _, t =
            a_slgmres_e(A, b; m_initial = 27, l = 3, d = 3, tol = 1e-12, maxit = 1000)
        @test flag
        @test relresvec[end] < 1e-12
        @test t > 0
    end

    @testset "Sherman5 grows m and converges" begin
        # Regression/lock test on the same problem and parameters used in
        # Cabral, Schaerer & Bhaya (2020)'s own numerical experiments
        # (master_algoritmos.m: mApd=28, dApd=2, lApd=2, alpha=2,
        # delta=0.8, epsilon=0.01). This port does not reproduce the
        # reference script's exact cycle count -- see the deviation
        # documented in a_slgmres_e.jl's docstring -- so this pins down
        # this port's own validated behavior (flag, final residual, and
        # that m does grow beyond m_initial) rather than the reference's
        # specific numbers.
        data_dir = joinpath(@__DIR__, "..", "..", "data")
        file = matopen(joinpath(data_dir, "sherman5.mat"))
        Problem = read(file, "Problem")
        close(file)
        A = Problem["A"]
        b = vec(Problem["b"])

        _, flag, relresvec, kdvec, t = a_slgmres_e(
            A,
            b;
            m_initial = 28,
            l = 2,
            d = 2,
            epsilon_threshold = 0.01,
            alpha_pd = [2.0, 0.8],
            tol = 1e-9,
            maxit = 1000,
        )
        @test flag
        @test relresvec[end] < 1e-9
        @test maximum(kdvec) > 28 + 2
        @test t > 0
    end

end

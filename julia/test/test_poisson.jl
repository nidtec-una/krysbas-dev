@testset "Poisson integration" begin

    @testset "1D Poisson, Dirichlet BC" begin
        f = x -> -6x
        g = x -> x^3
        NODES = 9
        INNERNODES = NODES - 2
        h = 1.0 / (NODES - 1)

        A = 2 * Matrix(I, INNERNODES, INNERNODES) -
            diagm(1 => ones(INNERNODES - 1)) -
            diagm(-1 => ones(INNERNODES - 1))
        xs = range(h, 1 - h; length=INNERNODES)
        b = h^2 * f.(xs)
        b[end] += g(1.0)
        uExact = g.(xs)

        u_pd, = pd_gmres(A, b; m_initial=3, tol=1e-9, maxit=100)
        u_lg, = lgmres(A, b; m=3, l=1, tol=1e-9, maxit=100)
        u_ge, = gmres_e(A, b; m=3, d=1, tol=1e-9, maxit=100)

        @test u_pd ≈ uExact atol=1e-6
        @test u_lg ≈ uExact atol=1e-6
        @test u_ge ≈ uExact atol=1e-6
    end

    @testset "1D Poisson, Dirichlet-Neumann BC" begin
        g = x -> (x + 1)^2
        NODES = 5
        INNERNODES = NODES - 2
        h = 1.0 / (NODES - 1)
        fx = -2.0
        DgN = 4.0

        n_sys = INNERNODES + 1
        A = 2 * Matrix(I, n_sys, n_sys) -
            diagm(1 => ones(n_sys - 1)) -
            diagm(-1 => ones(n_sys - 1))
        A[end, end] = 1.0

        b = h^2 * fx * ones(n_sys)
        b[1] += g(0.0)
        b[end] = 0.5 * h^2 * fx + h * DgN
        uExact = g.(range(h, 1.0; length=n_sys))

        u_pd, = pd_gmres(A, b; m_initial=3, tol=1e-9, maxit=100)
        u_lg, = lgmres(A, b; m=3, l=1, tol=1e-9, maxit=100)
        u_ge, = gmres_e(A, b; m=3, d=1, tol=1e-9, maxit=100)

        @test u_pd ≈ uExact atol=1e-6
        @test u_lg ≈ uExact atol=1e-6
        @test u_ge ≈ uExact atol=1e-6
    end

    @testset "1D Poisson, Robin-Robin BC" begin
        f = 2.0
        g = x -> (x + 1)^2
        NODES = 5
        alpha = 1.0; beta = -1.0; cStart = 3.0; cEnd = 0.0
        h = 1.0 / (NODES - 1)

        A = diagm(0 => fill(-2.0, NODES),
                  1 => ones(NODES - 1),
                  -1 => ones(NODES - 1))
        A[1, 1]     = h * alpha - 1
        A[end, end] = -h * beta - 1

        b = h^2 * f * ones(NODES)
        b[1]   = 0.5 * h^2 * f + h * cStart
        b[end] = 0.5 * h^2 * f - h * cEnd
        uExact = g.(range(0.0, 1.0; length=NODES))

        u_pd, = pd_gmres(A, b; m_initial=3, tol=1e-9, maxit=100)
        u_lg, = lgmres(A, b; m=3, l=1, tol=1e-9, maxit=100)
        # This 5×5 all-negative-eigenvalue matrix stagnates with SM harmonic
        # Ritz augmentation (m=3,d=1); d=0 dispatches to restarted GMRES(3).
        u_ge, = gmres_e(A, b; m=3, d=0, tol=1e-9, maxit=100)

        @test u_pd ≈ uExact atol=1e-6
        @test u_lg ≈ uExact atol=1e-6
        @test u_ge ≈ uExact atol=1e-6
    end

    @testset "2D Poisson, Dirichlet BC" begin
        NODES1D = 5
        INNODES1D = NODES1D - 2
        h = 1.0 / (NODES1D - 1)

        M = 4 * Matrix(I, INNODES1D, INNODES1D) -
            diagm(1 => ones(INNODES1D - 1)) -
            diagm(-1 => ones(INNODES1D - 1))
        N = -Matrix(I, INNODES1D, INNODES1D)
        Z = zeros(INNODES1D, INNODES1D)
        A = [M N Z; N M N; Z N M]

        b = zeros(INNODES1D * INNODES1D)
        h2 = h^2
        for i in 1:INNODES1D, j in 1:INNODES1D
            b[(i - 1) * INNODES1D + j] = h2 * 2 * pi^2 * sin(pi * i * h) * sin(pi * j * h)
        end

        uExact = [0.526514643772757; 0.744604150011472; 0.526514643772757;
                  0.744604150011472; 1.053029287545515; 0.744604150011472;
                  0.526514643772757; 0.744604150011472; 0.526514643772757]

        u_pd, = pd_gmres(A, b; m_initial=3, tol=1e-9, maxit=100)
        u_lg, = lgmres(A, b; m=3, l=1, tol=1e-9, maxit=100)
        u_ge, = gmres_e(A, b; m=3, d=1, tol=1e-9, maxit=100)

        @test u_pd ≈ uExact atol=1e-6
        @test u_lg ≈ uExact atol=1e-6
        @test u_ge ≈ uExact atol=1e-6
    end

    @testset "2D Poisson, left-to-right flow" begin
        M = [3.0 -1.0 0.0; -1.0 3.0 -1.0; 0.0 -1.0 3.0]
        N = -Matrix(I, 3, 3)
        P = [4.0 -1.0 0.0; -1.0 4.0 -1.0; 0.0 -1.0 4.0]
        Z = zeros(3, 3)
        A = [M N Z; N P N; Z N M]

        b = zeros(9)
        b[1] = 1.0; b[4] = 1.0; b[7] = 1.0

        uExact = [0.75; 0.50; 0.25; 0.75; 0.50; 0.25; 0.75; 0.50; 0.25]

        u_pd, = pd_gmres(A, b; m_initial=3, tol=1e-9, maxit=100)
        u_lg, = lgmres(A, b; m=3, l=1, tol=1e-9, maxit=100)
        u_ge, = gmres_e(A, b; m=3, d=1, tol=1e-9, maxit=100)

        @test u_pd ≈ uExact atol=1e-6
        @test u_lg ≈ uExact atol=1e-6
        @test u_ge ≈ uExact atol=1e-6
    end

end

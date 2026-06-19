@testset "pd_rule" begin

    # Common setup: n=100, mInitial=10, bounds [5, 50], step=2, αP=-3, αD=5
    n, m_init, m_min, m_max, m_step = 100, 10, 5, 50, 2
    αP, αD = -3.0, 5.0

    @testset "iter ≤ 2: returns m_initial unchanged" begin
        res = [1.0, 0.5]
        mj, _ = pd_rule(10, n, m_init, m_min, m_max, m_step, res, 2, αP, αD)
        @test mj == m_init
    end

    @testset "iter > 2: proportional term only" begin
        # mj = m + ceil(αP * (res[3]/res[2]))
        # αP = -3, res[3]/res[2] = 0.5/0.5 = 1.0 → ceil(-3.0) = -3
        # mj = 10 + (-3) = 7, which is ≥ m_min=5 → mj = 7
        res = [1.0, 0.5, 0.5]
        m = 10
        mj, _ = pd_rule(m, n, m_init, m_min, m_max, m_step, res, 3, αP, αD)
        expected = m + ceil(Int, αP * (res[3] / res[2]))
        @test mj == expected
    end

    @testset "iter > 3: full PD formula" begin
        res = [1.0, 0.5, 0.4, 0.2]
        m = 10
        mj, _ = pd_rule(m, n, m_init, m_min, m_max, m_step, res, 4, αP, αD)
        pd_term = αP * (res[4] / res[3]) + αD * ((res[4] - res[2]) / (2 * res[3]))
        expected = clamp(m + ceil(Int, pd_term), m_min, m_max)
        @test mj == expected
    end

    @testset "mj clamped to m_max" begin
        # Large positive residual ratio → mj would exceed m_max
        res = [1.0, 0.5, 0.4, 10.0]
        mj, _ = pd_rule(10, n, m_init, m_min, m_max, m_step, res, 4, 100.0, 0.0)
        @test mj == m_max
    end

    @testset "mj below m_min: m_initial steps up" begin
        # res[4]/res[3] = 0.5/0.4 = 1.25 → αP * 1.25 = -125 → ceil(-125) = -125
        # mj = 5 + (-125) = -120 < m_min=5 → step-up triggers
        res = [1.0, 0.5, 0.4, 0.5]
        m = 5
        mj, new_m_init = pd_rule(m, n, m_init, m_min, m_max, m_step, res, 4, -100.0, 0.0)
        @test new_m_init == m_init + m_step
        @test mj == m_init + m_step
    end

    @testset "output types are Int" begin
        res = [1.0, 0.5, 0.4, 0.2]
        mj, new_m_init = pd_rule(10, n, m_init, m_min, m_max, m_step, res, 4, αP, αD)
        @test mj isa Int
        @test new_m_init isa Int
    end

end

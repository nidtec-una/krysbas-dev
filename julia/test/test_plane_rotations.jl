@testset "plane_rotations" begin

    @testset "2×1 Hessenberg: known values" begin
        # H = [3; 4] (2×1), beta = 1.0
        # denom = 5, c = 3/5 = 0.6, s = 4/5 = 0.8
        # P*H = [5; 0],  P*g = [0.6; -0.8]
        H = reshape([3.0, 4.0], 2, 1)
        HUpTri, g = plane_rotations(H, 1.0)

        @test HUpTri[1, 1] ≈ 5.0
        @test abs(HUpTri[2, 1]) < 1e-14
        @test g[1] ≈ 0.6
        @test g[2] ≈ -0.8
    end

    @testset "subdiagonal of output is zero (upper triangular)" begin
        # 3×2 upper Hessenberg — after two Givens rotations every
        # subdiagonal entry must be (numerically) zero
        H = [2.0 1.0; 1.0 2.0; 0.0 1.0]
        HUpTri, _ = plane_rotations(H, 1.0)

        @test abs(HUpTri[2, 1]) < 1e-14
        @test abs(HUpTri[3, 2]) < 1e-14
    end

    @testset "Givens rotations preserve the norm of g" begin
        # Rotations are orthogonal: ‖g_out‖₂ = |beta|
        H = [2.0 1.0; 1.0 2.0; 0.0 1.0]
        beta = 3.7
        _, g = plane_rotations(H, beta)

        @test norm(g) ≈ beta
    end

    @testset "output dimensions match input" begin
        m = 4
        H = rand(m + 1, m)
        HUpTri, g = plane_rotations(H, 1.0)

        @test size(HUpTri) == (m + 1, m)
        @test length(g) == m + 1
    end

end

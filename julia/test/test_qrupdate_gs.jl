@testset "qrupdate_gs" begin

    @testset "full one-shot factorisation" begin
        rng = MersenneTwister(0)
        A = randn(rng, 10, 4)

        q_out, r_out = qrupdate_gs(A, zeros(0, 0), zeros(0, 0))

        @test q_out * r_out ≈ A
        @test q_out' * q_out ≈ I(4)
        @test all(abs.(tril(r_out, -1)) .< 1e-12)
    end

    @testset "R diagonal is non-negative" begin
        rng = MersenneTwister(1)
        A = randn(rng, 8, 3)

        _, r_out = qrupdate_gs(A, zeros(0, 0), zeros(0, 0))

        @test all(diag(r_out) .>= 0)
    end

    @testset "incremental matches full" begin
        rng = MersenneTwister(2)
        A = randn(rng, 10, 5)

        q_ref, r_ref = qrupdate_gs(A, zeros(0, 0), zeros(0, 0))

        q_inc = zeros(0, 0)
        r_inc = zeros(0, 0)
        for j = 1:5
            q_inc, r_inc = qrupdate_gs(A[:, 1:j], q_inc, r_inc)
        end

        @test q_inc ≈ q_ref
        @test r_inc ≈ r_ref
    end

    @testset "incremental orthogonality preserved" begin
        rng = MersenneTwister(3)
        A = randn(rng, 12, 6)

        q_inc = zeros(0, 0)
        r_inc = zeros(0, 0)
        for j = 1:6
            q_inc, r_inc = qrupdate_gs(A[:, 1:j], q_inc, r_inc)
            @test q_inc' * q_inc ≈ I(j)
        end
    end

    @testset "row growth between calls (the gmres_dr use case)" begin
        # In gmres_dr the Hessenberg gains one row at every Arnoldi step: at
        # step j the matrix is (j+1)-by-j. A must be upper Hessenberg
        # (A[i,j] == 0 for i > j+1), matching the actual gmres_dr use case:
        # qrupdate_gs extends the existing basis assuming the new bottom
        # row is zero for previously-seen columns, which only holds for
        # Hessenberg-structured input, not a fully dense matrix.
        rng = MersenneTwister(4)
        m = 8
        A = triu(randn(rng, m + 1, m), -1)

        q_inc = zeros(0, 0)
        r_inc = zeros(0, 0)
        for j = 1:m
            a_sub = A[1:(j+1), 1:j]
            q_inc, r_inc = qrupdate_gs(a_sub, q_inc, r_inc)

            @test q_inc * r_inc ≈ a_sub
            @test q_inc' * q_inc ≈ I(j)
        end
    end

    @testset "row growth residual estimate" begin
        rng = MersenneTwister(5)
        m = 6
        H = triu(randn(rng, m + 1, m), -1)
        c = randn(rng, m + 1)

        q_inc = zeros(0, 0)
        r_inc = zeros(0, 0)
        for j = 1:m
            h_sub = H[1:(j+1), 1:j]
            q_inc, r_inc = qrupdate_gs(h_sub, q_inc, r_inc)
        end

        d_inc = r_inc \ (q_inc' * c)
        res_inc = norm(c - H * d_inc)

        d_ref = pinv(H) * c
        res_ref = norm(c - H * d_ref)

        @test res_inc ≈ res_ref
    end

    @testset "single column" begin
        a = reshape([3.0, 4.0], 2, 1)

        q_out, r_out = qrupdate_gs(a, zeros(0, 0), zeros(0, 0))

        @test r_out[1, 1] ≈ 5.0
        @test q_out ≈ reshape([3.0, 4.0] / 5, 2, 1)
        @test (q_out'*q_out)[1, 1] ≈ 1.0
    end

end

@testset "modified_gram_schmidt_arnoldi" begin

    @testset "output dimensions" begin
        n, m = 10, 4
        A = diagm(0 => fill(2.0, n), 1 => fill(-1.0, n - 1), -1 => fill(-1.0, n - 1))
        v = ones(n) / sqrt(n)
        H, V, m_upd = modified_gram_schmidt_arnoldi(A, v, m)

        @test size(H) == (m + 1, m)
        @test size(V) == (n, m)
        @test m_upd == m
    end

    @testset "V columns are orthonormal" begin
        # Diagonal matrix with distinct eigenvalues avoids the near-breakdown
        # that arises with centrosymmetric matrices and uniform starting vectors.
        n, m = 8, 5
        A = diagm(collect(1.0:n))
        v = ones(n) / sqrt(n)
        _, V, _ = modified_gram_schmidt_arnoldi(A, v, m)

        @test V' * V ≈ I atol = 1e-12
    end

    @testset "happy breakdown: A proportional to I" begin
        # A = 2I with v = e₁ → Arnoldi breaks down at step 1:
        # W = 2e₁, H[1,1] = 2, W becomes 0, H[2,1] = 0 → early return.
        n = 4
        A = 2.0 * Matrix(I, n, n)
        v = [1.0, 0.0, 0.0, 0.0]
        H, V, m_upd = modified_gram_schmidt_arnoldi(A, v, n)

        @test m_upd == 1
        @test size(H) == (2, 1)
        @test size(V) == (n, 1)
        @test H[1, 1] ≈ 2.0
        @test abs(H[2, 1]) < 1e-15
    end

    @testset "Arnoldi residual relation" begin
        # At each step j, the MGS process guarantees:
        # norm(A*V[:,j] - sum_{i=1}^{j} H[i,j]*V[:,i]) == H[j+1, j]
        n, m = 6, 4
        A = diagm(collect(1.0:n))
        v = ones(n) / sqrt(n)
        H, V, _ = modified_gram_schmidt_arnoldi(A, v, m)

        for j = 1:m
            r = A * V[:, j]
            for i = 1:j
                r -= H[i, j] * V[:, i]
            end
            @test norm(r) ≈ H[j+1, j] atol = 1e-10
        end
    end

end

@testset "augmented_gram_schmidt_arnoldi" begin

    @testset "output dimensions" begin
        n, m, k = 10, 4, 2
        A = diagm(collect(1.0:n))
        v = ones(n) / sqrt(n)
        appendV = Matrix(I, n, k)  # first k standard basis vectors
        H, V, s = augmented_gram_schmidt_arnoldi(A, v, m, appendV)

        @test s == m + k
        @test size(H) == (s + 1, s)
        @test size(V) == (n, s)
    end

    @testset "V columns are orthonormal" begin
        n, m, k = 8, 3, 2
        A = diagm(collect(1.0:n))
        v = ones(n) / sqrt(n)
        appendV = Matrix(I, n, k)
        _, V, _ = augmented_gram_schmidt_arnoldi(A, v, m, appendV)

        @test V' * V ≈ I atol = 1e-12
    end

    @testset "k=0 matches standard MGS Arnoldi" begin
        n, m = 8, 4
        A = diagm(collect(1.0:n))
        v = ones(n) / sqrt(n)
        H1, V1, s1 = augmented_gram_schmidt_arnoldi(A, v, m, zeros(n, 0))
        H2, V2, s2 = modified_gram_schmidt_arnoldi(A, v, m)

        @test s1 == s2
        @test H1 ≈ H2 atol = 1e-14
        @test V1 ≈ V2 atol = 1e-14
    end

    @testset "happy breakdown with k=0" begin
        n, m = 4, 3
        A = 2.0 * Matrix(I, n, n)
        v = [1.0, 0.0, 0.0, 0.0]
        H, V, s = augmented_gram_schmidt_arnoldi(A, v, m, zeros(n, 0))

        @test s == 1
        @test size(V) == (n, 1)
        @test H[1, 1] ≈ 2.0
    end

    @testset "augmented columns increase subspace size" begin
        # With k extra vectors appended, s = m + k > m
        n, m, k = 10, 3, 3
        A = diagm(collect(1.0:n))
        v = ones(n) / sqrt(n)
        appendV = Matrix(I, n, k)
        _, _, s_aug = augmented_gram_schmidt_arnoldi(A, v, m, appendV)
        _, _, s_std = modified_gram_schmidt_arnoldi(A, v, m)

        @test s_aug == s_std + k
    end

end

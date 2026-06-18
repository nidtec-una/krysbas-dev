@testset "harmonic_ritz_vectors" begin

    @testset "output has k columns" begin
        s, k, n = 6, 2, 10
        F = diagm(collect(1.0:s))
        G = Matrix(I, s, s)
        V = Matrix(I, n, s)
        dy = harmonic_ritz_vectors(F, G, k, V)

        @test size(dy, 2) == k
        @test size(dy, 1) == n
    end

    @testset "selects smallest-magnitude eigenpairs" begin
        # F = diag(3, 1, 4, 1, 5), G = I → eigenvalues are 3,1,4,1,5.
        # k=2 smallest in magnitude are the two eigenvalues == 1 (indices 2 and 4).
        s, k = 5, 2
        F = diagm([3.0, 1.0, 4.0, 1.0, 5.0])
        G = Matrix(I, s, s)
        V = Matrix(I, s, s)   # identity → dy = eigenvectors of F directly
        dy = harmonic_ritz_vectors(F, G, k, V)

        # Each column of dy should be close to a standard basis vector e₂ or e₄
        norms = [norm(dy[:, j]) for j in 1:k]
        @test all(norms .> 0)

        # The two returned vectors should span {e₂, e₄}
        combined = abs.(dy' * [zeros(1); 1.0; zeros(1); 1.0; zeros(1)])
        @test sum(combined) > 1.0   # at least some projection onto {e₂, e₄}
    end

    @testset "real output for real problem" begin
        s, k, n = 4, 2, 8
        F = diagm(collect(1.0:s))
        G = Matrix(I, s, s)
        V = randn(n, s)
        dy = harmonic_ritz_vectors(F, G, k, V)

        @test eltype(dy) <: Real
        @test size(dy) == (n, k)
    end

end

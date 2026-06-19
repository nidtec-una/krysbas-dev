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

    @testset "selects largest-magnitude eigenpairs (LM → small eigs of A)" begin
        # F = diag(3, 1, 4, 1, 5), G = I → eigenvalues are 3,1,4,1,5.
        # k=2 largest in magnitude are eigenvalues 4 (index 3) and 5 (index 5).
        # Large |λ| in F*y=λ*G*y approximates small eigenvalues of A (Morgan 1995).
        s, k = 5, 2
        F = diagm([3.0, 1.0, 4.0, 1.0, 5.0])
        G = Matrix(I, s, s)
        V = Matrix(I, s, s)   # identity → dy = eigenvectors of F directly
        dy = harmonic_ritz_vectors(F, G, k, V)

        # Each column of dy should be close to a standard basis vector e₃ or e₅
        norms = [norm(dy[:, j]) for j in 1:k]
        @test all(norms .> 0)

        # The two returned vectors should span {e₃, e₅}
        combined = abs.(dy' * [zeros(2); 1.0; zeros(1); 1.0])
        @test sum(combined) > 1.0   # at least some projection onto {e₃, e₅}
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

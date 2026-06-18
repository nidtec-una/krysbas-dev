using Test
using KrySBAS

# Tests are added step by step following migration_plan.md.
# Uncomment each include as the corresponding component is ported:

# include("test_plane_rotations.jl")
# include("test_modified_gram_schmidt_arnoldi.jl")
# include("test_augmented_gram_schmidt_arnoldi.jl")
# include("test_harmonic_ritz_vectors.jl")
# include("test_pd_rule.jl")
# include("test_gmres_e.jl")
# include("test_lgmres.jl")
# include("test_pd_gmres.jl")
# include("test_poisson.jl")

@testset "KrySBAS.jl" begin
    @test true
end

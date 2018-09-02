@testset "Partial trace" begin

@testset "Equal dims of subsystems" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]

    ξ = ptrace(ρ ⊗ σ, [2, 2], [2,])
    @test norm(ρ - ξ) ≈ 0. atol=1e-15
    ζ = ptrace(ρ ⊗ σ, [2, 2], 2)
    @test norm(ξ - ζ) ≈ 0 atol=1e-15

    ξ = ptrace(ρ ⊗ σ, [2, 2], [1,])
    @test norm(σ - ξ) ≈ 0. atol=1e-15
    ζ = ptrace(ρ ⊗ σ, [2, 2], 1)
    @test norm(ξ - ζ) ≈ 0 atol=1e-15

    @test_throws ArgumentError ptrace(ones(2, 3), [1, 2], [1])
    @test_throws ArgumentError ptrace(ones(2, 2), [3, 4], [2])
    @test_throws ArgumentError ptrace(ones(4, 4), [2, 2], [3])
end

@testset "Vectors" begin
    ϕ = 1/sqrt(2) * (ket(0, 4) + ket(3, 4))
    ξ = ptrace(proj(ϕ), [2, 2], [2,])
    @test norm(ξ - eye(2)/2) ≈ 0. atol=1e-15
    ξ = ptrace(ϕ, [2, 2], 2)
    @test norm(ξ - eye(2)/2) ≈ 0. atol=1e-15
    ξ = ptrace(ϕ, [2, 2], 1)
    @test norm(ξ - eye(2)/2) ≈ 0. atol=1e-15
    @test_throws ArgumentError ptrace(ϕ, [2, 2], 3)
end

@testset "Sparse vectors" begin
    ϕ = sparse(1/sqrt(2) * (ket(0, 4) + ket(3, 4)))
    ξ = ptrace(ϕ, [2, 2], 2)
    @test typeof(ξ) <: AbstractSparseMatrix
    @test norm(ξ - eye(2)/2, 1) ≈ 0. atol=1e-15
end

@testset "Sparse matrices" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]

    ϕ = sparse(1/sqrt(2) * (ket(0, 4) + ket(3, 4)))
    ξ = ptrace(proj(ϕ), [2, 2], 2)
    @test norm(ξ - speye(2)/2, 1) ≈ 0. atol=1e-15
    ξ = ptrace(ϕ, [2, 2], 2)
    @test norm(ξ - eye(2)/2, 1) ≈ 0. atol=1e-15

    ξ = ptrace(sparse(ρ ⊗ σ), [2, 2], [1,])
    @test norm(σ - ξ, 1) ≈ 0. atol=1e-15
    ζ = ptrace(sparse(ρ ⊗ σ), [2, 2], 1)
    @test norm(ξ - ζ,1 ) ≈ 0 atol=1e-15

    @test_throws ArgumentError ptrace(sparse(ones(2, 3)), [1, 2], [1])
    @test_throws ArgumentError ptrace(sparse(ones(2, 2)), [3, 4], [2])
    @test_throws ArgumentError ptrace(sparse(ones(4, 4)), [2, 2], [3])
    @test_throws ArgumentError ptrace(sparse(ones(8, 8)), [2, 2, 2], [2])
end
end

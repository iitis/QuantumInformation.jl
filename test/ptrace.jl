@testset verbose=true "Partial trace" begin
using DoubleFloats: Double64, ComplexDF64

@testset verbose=true "Equal dims of subsystems" begin
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

@testset verbose=true "Vectors" begin
    ϕ = 1/sqrt(2) * (ket(1, 4) + ket(4, 4))
    ξ = ptrace(proj(ϕ), [2, 2], [2,])
    @test ξ ≈ I/2 atol=1e-15
    ξ = ptrace(ϕ, [2, 2], 2)
    @test ξ ≈ I/2 atol=1e-15
    ξ = ptrace(ϕ, [2, 2], 1)
    @test ξ ≈ I/2 atol=1e-15
    @test_throws ArgumentError ptrace(ϕ, [2, 2], 3)

    ψ = [sqrt(1 / 3), sqrt(2 / 3)]
    ϕ = [sqrt(2 / 3), sqrt(1 / 3)]
    @test ptrace(ψ⊗ϕ, [2, 2], 1) ≈ proj(ϕ) atol=1e-15
    @test ptrace(ψ⊗ϕ, [2, 2], 2) ≈ proj(ψ) atol=1e-15
end
@testset verbose=true "Double64 Support" begin
    ρ = Double64[0.25 0.25; 0.25 0.75]
    σ = Double64[0.4 0.1; 0.1 0.6] # Real matrix for simplicity, or complex
    
    # Test tensor product and ptrace
    big = kron(ρ, σ)
    # ptrace out second system -> should get ρ
    res = ptrace(big, [2, 2], 2)
    @test eltype(res) == Double64
    @test norm(res - ρ) ≈ 0.0 atol=1e-25
end
end

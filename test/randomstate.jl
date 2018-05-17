@testset "Random states" begin

@testset "random_ket" begin
    ϕ = zeros(Float64, 20)
    ψ = zeros(ComplexF64, 20)
    random_ket!(ϕ)
    random_ket!(ψ)
    @test length(ϕ) == 20
    @test length(ψ) == 20
    @test typeof(ϕ) == Vector{Float64}
    @test typeof(ψ) == Vector{ComplexF64}
    @test sum(abs2.(ϕ)) ≈ 1. atol=1e-15
    @test sum(abs2.(ψ)) ≈ 1. atol=1e-15

    ϕ = random_ket(Float64, 100)
    ψ = random_ket(ComplexF64, 100)
    @test length(ϕ) == 100
    @test length(ψ) == 100
    @test typeof(ϕ) == Vector{Float64}
    @test typeof(ψ) == Vector{ComplexF64}
    @test sum(abs2.(ϕ)) ≈ 1. atol=1e-15
    @test sum(abs2.(ψ)) ≈ 1. atol=1e-15

    ϕ = random_ket(100)
    @test length(ϕ) == 100
    @test typeof(ϕ) == Vector{ComplexF64}
    @test sum(abs2.(ϕ)) ≈ 1. atol=1e-15
end

@testset "random_mixed_state_hs" begin
    ρ = zeros(Float64, 20, 20)
    σ = zeros(ComplexF64, 20, 20)
    random_mixed_state!(ρ)
    random_mixed_state!(σ)
    @test size(ρ) == (20, 20)
    @test size(σ) == (20, 20)
    @test typeof(ρ) == Matrix{Float64}
    @test typeof(σ) == Matrix{ComplexF64}
    @test norm(ρ - ρ') ≈ 0. atol=1e-13 # is close to hermitian
    @test norm(σ - σ') ≈ 0. atol=1e-13 # is close to hermitian
    @test trace(ρ) ≈ 1. atol=1e-15
    @test trace(σ) ≈ 1. atol=1e-15

    ρ = random_mixed_state(Float64, 20)
    σ = random_mixed_state(ComplexF64, 20)
    @test size(ρ) == (20, 20)
    @test size(σ) == (20, 20)
    @test typeof(ρ) == Matrix{Float64}
    @test typeof(σ) == Matrix{ComplexF64}
    @test norm(ρ - ρ') ≈ 0. atol=1e-13 # is close to hermitian
    @test norm(σ - σ') ≈ 0. atol=1e-13 # is close to hermitian
    @test trace(ρ) ≈ 1. atol=1e-15
    @test trace(σ) ≈ 1. atol=1e-15

    ρ = random_mixed_state(20)
    @test size(ρ) == (20, 20)
    @test typeof(ρ) == Matrix{Complex{Float64}}
    @test norm(ρ - ρ') ≈ 0. atol=1e-13 # is close to hermitian
    @test trace(ρ) ≈ 1. atol=1e-15

    ρ = random_mixed_state(20, 0.5)
    @test size(ρ) == (20, 20)
    @test typeof(ρ) == Matrix{Complex{Float64}}
    @test norm(ρ - ρ') ≈ 0. atol=1e-13 # is close to hermitian
    @test trace(ρ) ≈ 1. atol=1e-15

    ρ = random_mixed_state(20, 10)
    @test size(ρ) == (20, 20)
    @test typeof(ρ) == Matrix{Complex{Float64}}
    @test norm(ρ - ρ') ≈ 0. atol=1e-13 # is close to hermitian
    @test trace(ρ) ≈ 1. atol=1e-15

    ρ = random_mixed_state(20, 0.3333333)
    @test size(ρ) == (20, 20)
    @test typeof(ρ) == Matrix{Complex{Float64}}
    @test norm(ρ - ρ') ≈ 0. atol=1e-13 # is close to hermitian
    @test trace(ρ) ≈ 1. atol=1e-15
end

@testset "random_jamiolkowski_state" begin
    n = 10
    J = zeros(Float64, n^2, n^2)
    random_jamiolkowski_state!(J)
    tr = ptrace(J, [n, n], [1])
    @test norm(tr - eye(n)/n) ≈ 0. atol=1e-13
    @test trace(J) ≈ 1. atol=1e-13
    @test typeof(J) == Matrix{Float64}
    J = zeros(ComplexF64, n^2, n^2)
    random_jamiolkowski_state!(J)
    tr = ptrace(J, [n, n], [1])
    @test norm(tr - eye(n)/n) ≈ 0. atol=1e-13
    @test trace(J) ≈ 1. atol=1e-13
    @test typeof(J) == Matrix{ComplexF64}
    J = random_jamiolkowski_state(Float64, n)
    tr = ptrace(J, [n, n], [1])
    @test size(J) == (n^2, n^2)
    @test trace(J) ≈ 1. atol=1e-13
    @test norm(tr - eye(n)/n) ≈ 0. atol=1e-14
    @test typeof(J) == Matrix{Float64}
    J = random_jamiolkowski_state(ComplexF64, n)
    tr = ptrace(J, [n, n], [1])
    @test size(J) == (n^2, n^2)
    @test trace(J) ≈ 1. atol=1e-13
    @test norm(tr - eye(n)/n) ≈ 0. atol=1e-14
    @test typeof(J) == Matrix{ComplexF64}
end

end

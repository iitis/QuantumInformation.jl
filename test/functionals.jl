@testset "Functiopnals" begin
@testset "trace_distance" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]

    @test trace_distance(ρ, ρ) ≈ 0 atol=1e-14
    @test_broken trace_distance(ρ, σ) ≈ 0.42426406871192857 atol=1e-15
end

@testset "fidelity_sqrt" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]
    @test fidelity_sqrt(ρ, σ) ≈ real(trace(sqrtm(sqrtm(ρ) * σ * sqrtm(ρ)))) atol=1e-15
end

@testset "fidelity" begin
    ϕ = ket(0, 2)
    ψ = ket(1, 2)
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]
    @test fidelity(ϕ, ϕ) ≈ 1. atol=1e-15
    @test fidelity(ϕ, ψ) ≈ 0. atol=1e-15
    @test fidelity(ρ, ψ) ≈ 0.75 atol=1e-15
    @test fidelity(ϕ, σ) ≈ 0.4 atol=1e-15
    @test fidelity(ρ, σ) ≈ real(trace(sqrtm(sqrtm(ρ) * σ * sqrtm(ρ))))^2 atol=1e-15
end

@testset "shannon_entropy" begin
    @test shannon_entropy(1/2) ≈ log(2) atol=1e-15
    @test shannon_entropy(1/4) ≈ 0.5623351446188083 atol=1e-15
    @test shannon_entropy(0.5*ones(20)) ≈ 10log(2) atol=1e-15
end

@testset "entropy" begin
    ϕ = ket(0, 2)
    ρ = [0.25 0.25im; -0.25im 0.75]
    @test entropy(ϕ) == 0
    #TODO: add tests for mixed states
end
end

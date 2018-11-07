@testset "Basic functions" begin

@testset "ket" begin
    ϕ = ket(1, 4)
    ψ = ComplexF64[1, 0, 0, 0]
    @test norm(ϕ - ψ) ≈ 0.

    @test typeof(ket(Float64, 1, 4)) == Vector{Float64}
    @test typeof(ket(ComplexF64, 1, 4)) == Vector{ComplexF64}

    @test_throws ArgumentError ket(4, 3)
end

@testset "bra" begin
    ϕ = bra(1, 4)
    ψ = ComplexF64[1 0 0 0]
    @test norm(ϕ - ψ) ≈ 0.

    @test_throws ArgumentError bra(4,3)

    @test typeof(bra(Float64, 1, 4)) == LinearAlgebra.Adjoint{Float64,Array{Float64,1}}
    @test typeof(bra(ComplexF64, 1, 4)) == LinearAlgebra.Adjoint{Complex{Float64},Array{Complex{Float64},1}}
end

@testset "ketbra" begin
    ϕψ = ketbra(1, 1, 4)
    αβ = zeros(ComplexF64, 4, 4)
    αβ[1, 1] = 1
    @test norm(ϕψ - αβ) ≈ 0.

    @test_throws ArgumentError ketbra(4, 4, 3)

    @test typeof(ketbra(Float64, 1, 1, 4)) == Matrix{Float64}
    @test typeof(ketbra(ComplexF64, 1, 1, 4)) == Matrix{ComplexF64}
end

@testset "proj" begin
    ϕ = ket(1, 4)
    ϕϕ = proj(ϕ)
    ψψ = zeros(ComplexF64, 4 ,4)
    ψψ[1, 1] = 1
    @test norm(ϕϕ - ψψ) ≈ 0.
end

@testset "res" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    ϕ = res(ρ)
    ψ = [0.25, 0.25im, -0.25im, 0.75]
    @test norm(ϕ - ψ) ≈ 0.
end

@testset "unres" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = unres(res(ρ))
    @test norm(ρ - σ) ≈ 0.

    expected = [
    1	2;
    3	4;
    5	6
    ]
    @test unres(1:6, 2) == expected
end


@testset "max_mixed" begin
    d = 10
    ρ = max_mixed(d)
    @test all(diag(ρ) .≈ 1/d)
end

@testset "max_entangled" begin
    ϕ = max_entangled(4)
    @test norm(ϕ) ≈ 1
    @test ϕ[1] ≈ 1/sqrt(2) atol=1e-15
    @test ϕ[4] ≈ 1/sqrt(2) atol=1e-15
    @test ϕ[2] ≈ 0 atol=1e-15
    @test ϕ[3] ≈ 0 atol=1e-15
end

@testset "werner_state" begin
    ρ = werner_state(4, 0.2222)
    @test tr(ρ) ≈ 1
    @test ishermitian(ρ)

    @test_throws ArgumentError werner_state(4, 1.2)
end
end

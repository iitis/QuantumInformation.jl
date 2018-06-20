@testset "Basic functions" begin

@testset "ket" begin
    ϕ = ket(0, 4)
    @test_throws ArgumentError ket(4,3)
    ψ = ComplexF64[1, 0, 0, 0]
    @test norm(ϕ - ψ) ≈ 0.
    # TODO : Fix these tests: types depend on julia version
    # @test typeof(ket(Float64, 0, 4)) == Vector{Float64}
    # @test typeof(ket(ComplexF64, 0, 4)) == Vector{ComplexF64}

    ϕ = ket(0, 4, sparse=true)
    @test ϕ[1] == 1
    @test size(ϕ) == (4,)
    @test nnz(ϕ) == 1

    @test_throws ArgumentError ket(4, 3, sparse=true)
end

@testset "bra" begin
    ϕ = bra(0, 4)
    ψ = ComplexF64[1 0 0 0]
    @test_throws ArgumentError bra(4,3)
    @test norm(ϕ - ψ) ≈ 0.

    ϕ = bra(0, 4, sparse=true)
    @test ϕ[1] == 1
    @test size(ϕ) == (1, 4)
    @test nnz(ϕ') == 1
    @test_throws ArgumentError ket(4, 3, sparse=true)
    # TODO : Fix these tests: types depend on julia version
    # @test typeof(bra(Float64, 0, 4)) == LinearAlgebra.Adjoint{Float64,Array{Float64,1}}
    # @test typeof(bra(ComplexF64, 0, 4)) == LinearAlgebra.Adjoint{Complex{Float64},Array{Complex{Float64},1}}
end

@testset "ketbra" begin
    ϕψ = ketbra(0, 0, 4)
    αβ = zeros(ComplexF64, 4, 4)
    αβ[1, 1] = 1
    @test norm(ϕψ - αβ) ≈ 0.
    @test_throws ArgumentError ketbra(4,4,3)

    ϕψ = ketbra(0, 0, 4, sparse=true)
    αβ = spzeros(ComplexF64, 4, 4)
    αβ[1, 1] = 1
    @test norm(ϕψ - αβ, 1) ≈ 0.
    @test size(ϕψ) == (4, 4)
    @test nnz(ϕψ) == 1
    @test_throws ArgumentError ketbra(4, 4, 3, sparse=true)
    # TODO : Fix these tests: types depend on julia version
    # @test typeof(ketbra(Float64, 0, 0, 4)) == Matrix{Float64}
    # @test typeof(ketbra(ComplexF64, 0, 0, 4)) == Matrix{ComplexF64}
end

@testset "proj" begin
    ϕ = ket(0, 4)
    ϕϕ = proj(ϕ)
    ψψ = zeros(ComplexF64, 4 ,4)
    ψψ[1, 1] = 1
    @test norm(ϕϕ - ψψ) ≈ 0.
end

@testset "base_matrices" begin
    d = 4
    m = collect(Matrix{ComplexF64}, base_matrices(4))
    for i=1:d, j=1:d
        v = trace(m[i]' * m[j])
        i == j ? @test(v == 1.) : @test(v == 0.)
    end
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
    a = [1 2.1 3; 4 5 6]
    @test norm(unres([1, 2.1, 3, 4 ,5, 6], 2, 3) - a) ≈ 0.

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

    ρ = max_mixed(d, sparse=true)
    @test typeof(ρ) <: AbstractSparseMatrix
    @test nnz(ρ) == d
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
    @test trace(ρ) ≈ 1
    @test ishermitian(ρ)

    @test_throws ArgumentError werner_state(4, 1.2)
end
end

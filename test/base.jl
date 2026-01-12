@testset verbose=true "Basic functions" begin
    using DoubleFloats: Double64, ComplexDF64
    using SparseArrays

    @testset verbose=true "ket" begin
        ϕ = ket(1, 4)
        ψ = ComplexF64[1, 0, 0, 0]
        @test norm(ϕ - ψ) ≈ 0.0

        @test typeof(ket(Vector{Float64}, 1, 4)) == Vector{Float64}
        @test typeof(ket(Vector{ComplexF64}, 1, 4)) == Vector{ComplexF64}

        @test typeof(ket(Float32, 1, 4)) == Vector{Float32}

        @test_throws ArgumentError ket(4, 3)
        @test_throws ArgumentError ket(-1, 2)
    end

    @testset verbose=true "bra" begin
        ϕ = bra(1, 4)
        ψ = ComplexF64[1 0 0 0]
        @test norm(ϕ - ψ) ≈ 0.0

        @test_throws ArgumentError bra(4, 3)

        @test typeof(bra(Float32, 1, 4)) ==
              LinearAlgebra.Adjoint{Float32, Array{Float32, 1}}

        @test typeof(bra(Vector{Float64}, 1, 4)) ==
              LinearAlgebra.Adjoint{Float64, Array{Float64, 1}}
        @test typeof(bra(Vector{ComplexF64}, 1, 4)) ==
              LinearAlgebra.Adjoint{Complex{Float64}, Array{Complex{Float64}, 1}}
    end

    @testset verbose=true "ketbra" begin
        ϕψ = ketbra(1, 1, 4)
        αβ = zeros(ComplexF64, 4, 4)
        αβ[1, 1] = 1
        @test norm(ϕψ - αβ) ≈ 0.0

        ψ = ketbra(3, 2, 3, 4) * ket(2, 3)
        @test ψ[3] ≈ 1

        @test_throws ArgumentError ketbra(4, 4, 3)
        @test_throws ArgumentError ketbra(1, 1, -1)
        @test_throws ArgumentError ketbra(3, 1, 2)

        @test typeof(ketbra(Float32, 1, 1, 4)) == Matrix{Float32}

        @test typeof(ketbra(Matrix{Float64}, 1, 1, 4)) == Matrix{Float64}
        @test typeof(ketbra(Matrix{ComplexF64}, 1, 1, 4)) == Matrix{ComplexF64}
    end

    @testset verbose=true "proj" begin
        ϕ = ket(1, 4)
        ϕϕ = proj(ϕ)
        ψψ = zeros(ComplexF64, 4, 4)
        ψψ[1, 1] = 1
        @test norm(ϕϕ - ψψ) ≈ 0.0
    end

    @testset verbose=true "bloch_vector" begin
        @testset verbose=true "Test 1: Identity matrix" begin
            ρ = Matrix{ComplexF64}(I, 2, 2)
            @test bloch_vector(ρ) ≈ ComplexF64[0.0, 0.0, 1.0]
        end

        @testset verbose=true "Test 2: X basis state" begin
            ρ = [0.5 0.5; 0.5 0.5]
            @test bloch_vector(ρ) ≈ Float64[1.0, 0.0, 0.0]
        end

        @testset verbose=true "Test 3: Y basis state" begin
            ρ = [0.5 0.5im; -0.5im 0.5]
            @test bloch_vector(ρ) ≈ ComplexF64[0.0, 1.0, 0.0]
        end

        @testset verbose=true "Test 4: Z basis state" begin
            ρ = [1.0 0.0; 0.0 0.0]
            @test bloch_vector(ρ) ≈ Float64[0.0, 0.0, 1.0]
        end

        @testset verbose=true "Test 5: Arbitrary density matrix" begin
            ρ = [0.6 0.3-0.2im; 0.3+0.2im 0.4]
            @test bloch_vector(ρ) ≈ ComplexF64[0.6, -0.4, 0.2]
        end
    end

    @testset verbose=true "res" begin
        ρ = [0.25 0.25im; -0.25im 0.75]
        ϕ = res(ρ)
        ψ = [0.25, 0.25im, -0.25im, 0.75]
        @test norm(ϕ - ψ) ≈ 0.0
    end

    @testset verbose=true "unres" begin
        ρ = [0.25 0.25im; -0.25im 0.75]
        σ = unres(res(ρ))
        @test norm(ρ - σ) ≈ 0.0

        expected = [
            1 2;
            3 4;
            5 6
        ]
        @test unres(1:6, 2) == expected
    end

    @testset verbose=true "max_mixed" begin
        d = 10
        ρ = max_mixed(d)
        @test all(diag(ρ) .≈ 1/d)
    end

    @testset verbose=true "max_entangled" begin
        ϕ = max_entangled(4)
        @test norm(ϕ) ≈ 1
        @test ϕ[1] ≈ 1/sqrt(2) atol=1e-15
        @test ϕ[4] ≈ 1/sqrt(2) atol=1e-15
        @test ϕ[2] ≈ 0 atol=1e-15
        @test ϕ[3] ≈ 0 atol=1e-15
    end

    @testset verbose=true "werner_state" begin
        ρ = werner_state(4, 0.2222)
        @test tr(ρ) ≈ 1
        @test ishermitian(ρ)

        @test_throws ArgumentError werner_state(4, 1.2)
        @test_throws ArgumentError werner_state(4, -0.1)
    end
    @testset verbose=true "Double64 Support" begin
        @testset verbose=true "ket" begin
            ϕ = ket(ComplexDF64, 1, 4)
            @test eltype(ϕ) == ComplexDF64
            @test norm(ϕ) ≈ 1.0
        end

        @testset verbose=true "bra" begin
            ϕ = bra(ComplexDF64, 1, 4)
            @test eltype(ϕ) == ComplexDF64
            @test norm(ϕ) ≈ 1.0
        end

        @testset verbose=true "max_entangled" begin
            ϕ = max_entangled(ComplexDF64, 4)
            @test eltype(ϕ) == ComplexDF64
            @test norm(ϕ) ≈ 1.0
        end

        @testset verbose=true "werner_state" begin
            ρ = werner_state(4, Double64(0.2))
            @test eltype(ρ) == ComplexDF64
            @test tr(ρ) ≈ 1.0
        end
    end

    @testset verbose=true "Sparse Support" begin
        @testset verbose=true "ket" begin
            ϕ = ket(SparseVector{ComplexF64, Int}, 1, 4)
            @test ϕ isa SparseVector
            @test ϕ[1] == 1.0 + 0.0im
        end

        @testset verbose=true "bra" begin
            ϕ = bra(SparseVector{ComplexF64, Int}, 1, 4)
            @test ϕ' isa SparseVector # bra returns Adjoint, transpose is the vector
        end

        @testset verbose=true "ketbra" begin
            ρ = ketbra(SparseMatrixCSC{ComplexF64, Int}, 1, 1, 2)
            @test ρ isa SparseMatrixCSC
            @test ρ[1, 1] == 1.0 + 0.0im
        end

        @testset verbose=true "max_entangled" begin
            ϕ = max_entangled(SparseVector{ComplexF64, Int}, 4)
            @test ϕ isa SparseVector
            @test norm(ϕ) ≈ 1.0
        end
    end
end

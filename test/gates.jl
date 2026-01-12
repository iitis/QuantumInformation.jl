@testset verbose=true "Gates" begin
    using DoubleFloats: Double64, ComplexDF64
    using SparseArrays

    @testset verbose=true "QFT" begin
        d = 10
        F = qft(d)
        @test size(F) == (d, d)
        @test F'*F ≈ I atol=1e-13
        @test norm(abs.(F) - fill(1/sqrt(d), d, d)) ≈ 0 atol=1e-13
    end

    @testset verbose=true "Grover" begin
        d = 10
        G = grover(d)
        @test G'*G ≈ I atol=1e-13
        @test size(G) == (d, d)
    end

    @testset verbose=true "hadamard" begin
        d = 16
        H = hadamard(d)
        @test size(H) == (d, d)
        @test H'*H ≈ I atol=1e-13
        @test H' ≈ H atol=1e-13
        @test_throws ArgumentError hadamard(10)
    end

    @testset verbose=true "Pauli matrices" begin
        @test size(sx) == (2, 2)
        @test size(sy) == (2, 2)
        @test size(sz) == (2, 2)

        @test sx*sy - sy*sx == -2im * sz
    end

    @testset verbose=true "Double64 Support" begin
        d = 10

        @testset verbose=true "QFT" begin
            F = qft(ComplexDF64, d)
            @test eltype(F) == ComplexDF64
            @test F'*F ≈ I atol=1e-25
        end

        @testset verbose=true "Grover" begin
            G = grover(ComplexDF64, d)
            @test eltype(G) == ComplexDF64
            @test G'*G ≈ I atol=1e-25
        end

        @testset verbose=true "Hadamard" begin
            H = hadamard(Double64, 16)
            @test eltype(H) == Double64 # Note: hadamard returns real matrix
            @test H'*H ≈ I atol=1e-25
        end
    end
end

@testset verbose=true "Sparse Support" begin
    d = 4
    @testset verbose=true "Identity" begin
        id = 𝕀(SparseMatrixCSC{ComplexF64, Int}, d)
        @test id isa SparseMatrixCSC
        @test isapprox(id, I(d))
    end

    @testset verbose=true "QFT" begin
        F = qft(SparseMatrixCSC{ComplexF64, Int}, d)
        @test F isa SparseMatrixCSC
        @test F'*F ≈ I atol=1e-13
    end

    @testset verbose=true "Grover" begin
        G = grover(SparseMatrixCSC{ComplexF64, Int}, d)
        @test G isa SparseMatrixCSC
        @test G'*G ≈ I atol=1e-13
    end

    @testset verbose=true "Hadamard" begin
        H = hadamard(SparseMatrixCSC{Float64, Int}, 4)
        @test H isa SparseMatrixCSC
        @test H'*H ≈ I atol=1e-13
    end
end

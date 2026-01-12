@testset verbose=true "Utility functions" begin
using DoubleFloats: Double64, ComplexDF64
using SparseArrays

@testset verbose=true "number2mixedradix" begin
   number = 486
   bases = Int64[8, 42, 2]
   @test number2mixedradix(number, bases) == Int64[5, 33, 0]
end

@testset verbose=true "mixedradix2number" begin
   number = Int64[5, 33, 0]
   bases = Int64[8, 42, 2]
   @test mixedradix2number(number, bases) == 486
end

@testset verbose=true "renormalize" begin
    rng = MersenneTwister(1234);
    v = randn(rng, 10)
    renormalize!(v)

    @test norm(v) ≈ 1 atol=1e-13

    A = randn(rng, 10, 10)
    renormalize!(A)

    @test tr(A) ≈ 1 atol=1e-13
end

@testset verbose=true "funcmh" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    R = zero(ρ)
    ref = exp(ρ)
    funcmh!(x->exp(x), Hermitian(ρ), R)
    @test R ≈ ref atol=1e-15
    ρ = [0.25 0.25im; -0.25im 0.75]
    R = funcmh(x->exp(x), ρ)
    @test R ≈ ref atol=1e-15

    ρ = [0.25 0.25im; -0.25im 0.75]
    R = zero(ρ)
    funcmh!(x->exp(x), ρ, R)

    @test R ≈ ref atol=1e-15
end

    @testset verbose=true "Double64" begin
        v = Double64[1, 1]
        renormalize!(v)
        @test norm(v) ≈ 1.0 atol=1e-25
        @test v[1] isa Double64
        
        m = Double64[1 0; 0 1]
        renormalize!(m)
        @test tr(m) ≈ 1.0 atol=1e-25
        @test m[1,1] isa Double64
    end
    
    @testset verbose=true "Sparse Support" begin
        v = sparsevec([1], [2.0], 4)
        renormalize!(v)
        @test norm(v) ≈ 1.0
        @test v[1] ≈ 1.0
        
        m = sparse([1], [1], [2.0], 4, 4)
        renormalize!(m)
        @test tr(m) ≈ 1.0
        @test m[1,1] ≈ 1.0
        
        # Identity check
        @test isidentity(sparse(I, 4, 4))
        
        # Positive check (might covert to dense)
        @test ispositive(sparse(I, 4, 4))
    end
end

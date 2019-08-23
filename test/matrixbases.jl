@testset "MatrixBases" begin

@testset "HermitianBasisIterator" begin
    d = 4
    m = collect(HermitianBasisIterator{Matrix{ComplexF64}}(d))
    @test [tr(m[i]' * m[j]) for i=1:d, j=1:d] ≈ Matrix{Float64}(I, d, d)
end

@testset "ChannelBasisIterator" begin
    d1 = 2
    d2 = 2
    d = d1^2 * d2^2 - d1^2 + 1
    m = collect(ChannelBasisIterator{Matrix{ComplexF64}}(d1,d2))
    @test [tr(m[i]' * m[j]) for i=1:d, j=1:d] ≈ Matrix{Float64}(I, d, d)
end

@testset "represent, combine" begin
    d = 4
    A = reshape(collect(1:16), d, d) + reshape(collect(1:16), d, d)'
    vA = represent(HermitianBasis{Matrix{ComplexF64}}(d), A)
    Ap = combine(HermitianBasis{Matrix{ComplexF64}}(d), vA)
    @test A ≈ Ap
    B = A*A'
    vB = represent(HermitianBasis{Matrix{ComplexF64}}(d), B)
    Bp = combine(HermitianBasis{Matrix{ComplexF64}}(d), vB)
    @test B ≈ Bp

    vB = represent(HermitianBasis{Matrix{ComplexF32}}(d), B)
    @test eltype(vB) == Float32

    C = Float16[1 2; 3 4]
    C += C'
    vC = represent(HermitianBasis, C)
    @test eltype(vC) == eltype(C)
    @test length(vC) == prod(size(C))
end

@testset "represent, combine" begin
    d1 = 2
    d2 = 4
    A = reshape(collect(1:16), d1 * d2, d1) * reshape(collect(1:16), d1 * d2, d1)'
    B = Matrix{Float64}(I, d2, d2) ⊗ (ptrace(A, [d2, d1], 1))^(-1/2)
    A = B * A * B'
    vA = represent(ChannelBasis{Matrix{ComplexF64}}(d1, d2), A)
    Ap = combine(ChannelBasis{Matrix{ComplexF64}}(d1, d2), vA)
    @test A ≈ Ap.matrix
    
    A = reshape(collect(1:64), d1 * d2, d1 * d2) * reshape(collect(1:64), d1 * d2, d1 * d2)' + Matrix{Float64}(I, d1 * d2, d1 * d2)
    B = Matrix{Float64}(I, d1, d1) ⊗ (ptrace(A, [d1, d2], 1))^(-1/2)
    B = B * A * B'
    vB = represent(ChannelBasis{Matrix{ComplexF64}}(d2, d1), B)
    Bp = combine(ChannelBasis{Matrix{ComplexF64}}(d2, d1), vB)
    @test B ≈ Bp.matrix

    #vB = represent(ChannelBasis{Matrix{ComplexF32}}(d1,d2), B)
    #@test eltype(vB) == Float32

    # C = Float16[1 2; 3 4]
    # C += C'
    # vC = represent(ChannelBasis, C)
    # @test eltype(vC) == eltype(C)
    # @test length(vC) == prod(size(C))
end

@testset "hermitainbasis" begin
    @test hermitianbasis(Matrix{Float32}, 2) == HermitianBasisIterator{Matrix{Float32}}(2)
    @test hermitianbasis(2) == HermitianBasisIterator{Matrix{ComplexF64}}(2)
end

# @testset "channelbasis" begin
#    @test channelbasis(Matrix{Float32}, 2,2) == ChannelBasisIterator{Matrix{Float32}}(2,2)
#    @test channelbasis(2,2) == ChannelBasisIterator{Matrix{ComplexF64}}(2,2)
# end

end
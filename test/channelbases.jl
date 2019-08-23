@testset "MatrixBases" begin

@testset "ChannelBasisIterator" begin
    d1=2
    d2=2
    d=d1^2*d^2-d1^2+1
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

@testset "hermitainbasis" begin
    @test hermitianbasis(Matrix{Float32}, 2) == HermitianBasisIterator{Matrix{Float32}}(2)
    @test hermitianbasis(2) == HermitianBasisIterator{Matrix{ComplexF64}}(2)
end

end
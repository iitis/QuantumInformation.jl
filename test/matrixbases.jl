@testset "MatrixBases" begin

@testset "HermitianBasisIterator" begin
    d = 4
    m = collect(HermitianBasisIterator{Matrix{ComplexF64}}(d))
    @test [tr(m[i]' * m[j]) for i=1:d, j=1:d] ≈ Matrix{Float64}(I, d, d)
end

@testset "ChannelBasisIteratorsquare" begin
    idim = 2
    odim = idim 
    d = idim^2 * odim^2 - idim^2 + 1
    m = collect(ChannelBasisIterator{Matrix{ComplexF64}}(idim,odim))
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

@testset "representchannel, combinechannel" begin
    (idim, odim) = (2,4)
    A = reshape(collect(1:16), idim * odim, idim) * reshape(collect(1:16), idim * odim, idim)'
    B = Matrix{Float64}(I, odim, odim) ⊗ (ptrace(A, [odim, idim], 1))^(-1/2)
    A = B * A * B'
    vA = represent(channelbasis(Matrix{ComplexF64}, idim, odim), A)
    Ap = combine(channelbasis(Matrix{ComplexF64}, idim, odim), vA)
    @test A ≈ Ap.matrix
    
    A = reshape(collect(1:64), idim * odim, idim * odim) * reshape(collect(1:64), idim * odim, idim * odim)' + Matrix{Float64}(I, idim * odim, idim * odim)
    B = Matrix{Float64}(I, idim, idim) ⊗ (ptrace(A, [idim, odim], 1))^(-1/2)
    B = B * A * B'
    vB = represent(channelbasis(Matrix{ComplexF64}, odim, idim), B)
    Bp = combine(channelbasis(Matrix{ComplexF64}, odim, idim), vB)
    @test B ≈ Bp.matrix

    vB = represent(channelbasis(Matrix{ComplexF32}, idim,odim), B)
    @test eltype(vB) == Float64

end

@testset "hermitainbasis" begin
    @test hermitianbasis(Matrix{Float32}, 2) == HermitianBasisIterator{Matrix{Float32}}(2)
    @test hermitianbasis(2) == HermitianBasisIterator{Matrix{ComplexF64}}(2)
end

@testset "channelbasis" begin
    @test channelbasis(Matrix{Float32}, 2) == ChannelBasis{Matrix{Float32}}(2,2)
    @test channelbasis(2,3) ==  ChannelBasis{Matrix{ComplexF64}}(2,3)
end

end
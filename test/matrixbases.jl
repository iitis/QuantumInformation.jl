@testset verbose=true "MatrixBases" begin
using DoubleFloats: Double64, ComplexDF64
using SparseArrays

@testset verbose=true "HermitianBasisIterator" begin
    d = 4
    m = collect(HermitianBasisIterator{Matrix{ComplexF64}}(d))
    @test [tr(m[i]' * m[j]) for i=1:d, j=1:d] ≈ Matrix{Float64}(I, d, d)
end

@testset verbose=true "ChannelBasisIteratorsquare" begin
    idim = 2
    odim = idim 
    d = idim^2 * odim^2 - idim^2 + 1
    m = collect(ChannelBasisIterator{Matrix{ComplexF64}}(idim,odim))
    @test [tr(m[i]' * m[j]) for i=1:d, j=1:d] ≈ Matrix{Float64}(I, d, d)
end

@testset verbose=true "represent, combine" begin
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

@testset verbose=true "representchannel, combinechannel" begin
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

@testset verbose=true "hermitainbasis" begin
    @test hermitianbasis(Matrix{Float32}, 2) == HermitianBasisIterator{Matrix{Float32}}(2)
    @test hermitianbasis(2) == HermitianBasisIterator{Matrix{ComplexF64}}(2)
end

@testset verbose=true "channelbasis" begin
    @test channelbasis(Matrix{Float32}, 2) == ChannelBasis{Matrix{Float32}}(2,2)
    @test channelbasis(2,3) ==  ChannelBasis{Matrix{ComplexF64}}(2,3)
end

@testset verbose=true "Double64 Support" begin
    d = 2
    m = collect(HermitianBasisIterator{Matrix{ComplexDF64}}(d))
    @test eltype(m[1]) == ComplexDF64
    @test [tr(m[i]' * m[j]) for i=1:d, j=1:d] ≈ Matrix{Double64}(I, d, d) atol=1e-15

    v = represent(HermitianBasis{Matrix{ComplexDF64}}(d), m[1])
    @test eltype(v) == Double64
    combined = combine(HermitianBasis{Matrix{ComplexDF64}}(d), v)
    @test combined ≈ m[1] atol=1e-25
    v = represent(HermitianBasis{Matrix{ComplexDF64}}(d), m[1])
    @test eltype(v) == Double64
    combined = combine(HermitianBasis{Matrix{ComplexDF64}}(d), v)
    @test combined ≈ m[1] atol=1e-25
end

@testset verbose=true "Sparse Support" begin
    d = 2
    # Check HermitianBasis with SparseMatrixCSC
    basis = HermitianBasis{SparseMatrixCSC{ComplexF64, Int}}(d)
    m = collect(basis.iterator)
    @test m[1] isa SparseMatrixCSC
    @test [real(tr(m[i]' * m[j])) for i=1:d^2, j=1:d^2] ≈ Matrix{Float64}(I, d^2, d^2)

    # Represent sparse matrix
    A = sparse([0. 1.; 1. 0.])
    vA = represent(basis, A)
    Ap = combine(basis, vA)
    @test Ap isa SparseMatrixCSC
    @test A ≈ Ap

    # ChannelBasis
    idim, odim = 2, 2
    cbasis = channelbasis(SparseMatrixCSC{ComplexF64, Int}, idim, odim)
    m_chan = collect(cbasis.iterator)
    @test m_chan[1] isa SparseMatrixCSC
    
    # Check orthonormality
    n = length(cbasis.iterator)
    # This might be slow if n is large, d=2 -> n = 4*4 - 4 + 1 = 13. manageable.
    @test [real(tr(m_chan[i]' * m_chan[j])) for i=1:n, j=1:n] ≈ Matrix{Float64}(I, n, n)
    
    # Represent sparse channel (Process Matrix)
    # A simple channel: Identity channel. J = sum |ii><jj| ⊗ |i><j|?
    # Or just random sparse
    S = sparse(Matrix(I, idim*odim, idim*odim))
    vS = represent(cbasis, S)
    Sp = combine(cbasis, vS)
    
    @test Sp.matrix isa SparseMatrixCSC
    @test Sp.matrix ≈ S # limit_type removes nearly zero values
end

end
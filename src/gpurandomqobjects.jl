using CuArrays
using GPUArrays
import .CuRandomMatrices: curand

renormalize!(x::CuVector) = (x = x ./ CuArrays.norm(x))
renormalize!(x::CuMatrix) = (x = x ./ CuArrays.tr(x))
# this works only for positive matrices!!
# invsqrt
function invsqrt!(x::CuMatrix)
    F = CuArrays.svd!(x)
    S = 1 ./ sqrt.(F.S)
    (CuArrays.transpose(S) .* F.U) * F.U'
end

function invsqrt(x::CuMatrix)
    y = copy(x)
    invsqrt!(y)
end

# kron for CuArrays
function kron(a::CuMatrix{T1}, b::CuMatrix{T2}) where {T1<:Number, T2<:Number}
    dims = map(prod, zip(size(a), size(b)))
    T = promote_type(T1, T2)
    c = CuMatrix{T}(undef, dims...)
    gpu_call(c, (c, a, b)) do state, c, a, b
        (i, j) = @cartesianidx c
        (p, q) = size(b)
        k = (i - 1) ÷ p + 1
        l = (j - 1) ÷ q + 1
        m = i - ((i - 1) ÷ p) * p
        n = j - ((j - 1) ÷ q) * q
        @inbounds c[i, j] = a[k, l] * b[m, n]
        return
    end
    return c
end

function curand(h::HaarKet{1})
    ψ = CuArrays.randn(h.d)
    renormalize!(ψ)
 end

 function curand(h::HaarKet{2})
     ψ = CuArrays.randn(h.d) + 1im * CuArrays.randn(h.d)
     renormalize!(ψ)
  end

function curand(hs::HilbertSchmidtStates)
    ρ = curand(hs.w)
    renormalize!(ρ)
end

function curand(c::ChoiJamiolkowskiMatrices)
    z = curand(c.w)
    y = ptrace(z, [c.odim, c.idim], [1])
    sy = invsqrt(y)
    # onesy = CuMatrix{eltype(sy)}}(undef, size(sy) .* c.odim)
    # gpu_call(onesy, (onesy, sy)) do state, onesy, sy
    #     (i, j) = @cartesianidx onesy
    #     onesy[i, j] = ? 
    # end
    #TODO: fix this, can be accelerated with a custom kernel
    onesy = CuMatrix(I(c.odim)) ⊗ sy
    DynamicalMatrix(onesy * z * onesy, c.idim, c.odim)
end


function curand(c::HaarPOVM{N}) where N
    # TODO: this should be on the gpu in one go
    # TODO: use slicing??
    V = curand(c.c)
    POVMMeasurement([V'*(ketbra(CuMatrix{ComplexF32}, i, i, c.odim) ⊗ 𝕀(N))*V for i=1:c.odim])
end

function curand(c::VonNeumannPOVM)
    # TODO: this should be on the gpu in one go
    V = curand(c.c)
    POVMMeasurement([proj(V[:, i]) for i=1:c.d])
end

function curand(c::WishartPOVM)
    # TODO: this should be on the gpu in one go
    Ws = map(x->curand(x), c.c)
    S = sum(Ws)
    Ssq = invsqrt(S)
    POVMMeasurement([Ssq * W * Ssq for W=Ws])
end
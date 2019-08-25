using CuArrays
import .CuRandomMatrices: curand

renormalize!(x::CuVector) = (x = x ./ CuArrays.norm(x))
renormalize!(x::CuMatrix) = (x = x ./ CuArrays.tr(x))
# this works only for positive matrices!!
function invsqrt!(x::CuMatrix)
    F = CuArrays.svd!(x)
    S = 1 ./ sqrt.(F.S)
    (CuArrays.transpose(S) .* F.U) * F.U'
end

function invsqrt(x::CuMatrix)
    y = copy(x)
    invsqrt!(y)
end

function curand(h::HaarKet{1})
    ψ = CuArrays.randn(h.d)
    renormalize!(ψ)
 end

 function curand(h::HaarKet{2})
     ψ = CuArrays.randn(h.d) + 1im * CuArrays.randn(h.d)
     renormalize!(ψ)
  end

function curand(hs::HilbertSchmidtStates{β, K}) where {β, K}
    ρ = curand(hs.w)
    renormalize!(ρ)
end

function curand(c::ChoiJamiolkowskiMatrices{β, K}) where {β, K}
    error("Not iplmeneted")
    z = curand(c.w)
    y = ptrace(z, [c.odim, c.idim], [1])
    sy = invsqrt(y)
    onesy = Matrix(I, c.odim, c.odim) ⊗ sy # onesy = eye(c.odim) ⊗ sy
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
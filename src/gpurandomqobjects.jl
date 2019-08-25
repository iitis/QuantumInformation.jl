using CuArrays
using ..CuRandomMatrices

renormalize!(x::CuVector) = (x = x ./ CuArrays.norm(x))
renormalize!(x::CuMatrix) = (x = x ./ CuArrays.tr(x))

function curand(h::HaarKet{1})
    ψ = CuArrays.randn(h.d)
    renormalize!(ψ)
 end

 function curand(h::HaarKet{2})
     ψ = CuArrays.randn(h.d) + 1im * CuArrays.randn(h.d)
     renormalize!(ψ)
  end

function curand(hs::HilbertSchmidtStates{β, K}) where {β, K}
    ρ = CuRandomMatrices.curand(hs.w)
    renormalize!(ρ)
end

function curand(c::ChoiJamiolkowskiMatrices{β, K}) where {β, K}
    error("Not iplmeneted")
    z = curand(c.w)
    y = ptrace(z, [c.odim, c.idim], [1])
    sy = funcmh!(x -> 1 / sqrt(x), y)
    onesy = Matrix(I, c.odim, c.odim) ⊗ sy # onesy = eye(c.odim) ⊗ sy
    DynamicalMatrix(onesy * z * onesy, c.idim, c.odim)
end


function curand(c::HaarPOVM{N}) where N
    error("Not iplmeneted")
    V = curand(c.c)
    POVMMeasurement([V'*(ketbra(i, i, c.odim) ⊗ 𝕀(N))*V for i=1:c.odim])
end

function curand(c::VonNeumannPOVM)
    error("Not iplmeneted")
    V = curand(rng, c.c)
    POVMMeasurement([proj(V[:, i]) for i=1:c.d])
end

function curand(c::WishartPOVM)
    error("Not iplmeneted")
    Ws = map(x->curand(x), c.c)
    S = sum(Ws)
    Ssq = funcmh!(x->1/sqrt(x), S)
    POVMMeasurement([Ssq * W * Ssq for W=Ws])
end
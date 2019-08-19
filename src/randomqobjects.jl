export HaarKet, HilbertSchmidtStates, ChoiJamiolkowskiMatrices,
    HaarPOVM, WishartPOVM, VonNeumannPOVM


struct HaarKet{β} <: QIContinuousMatrixDistribution
    d::Int
end

HaarKet(d::Int) = HaarKet{2}(d)

function rand(rng::AbstractRNG, h::HaarKet{1})
    ψ = randn(rng, h.d)
    renormalize!(ψ)
    ψ
 end

 function rand(rng::AbstractRNG, h::HaarKet{2})
     ψ = randn(rng, h.d) + 1im * randn(rng, h.d)
     renormalize!(ψ)
     ψ
  end

# Random mixed states
struct HilbertSchmidtStates{β, K} <: QIContinuousMatrixDistribution
    w::WishartEnsemble
    d::Int

    function HilbertSchmidtStates{β, K}(d::Int) where {β, K}
        w = WishartEnsemble{β, K}(d)
        new(w, w.d)
    end
end
HilbertSchmidtStates{β}(d::Int) where β = HilbertSchmidtStates{β, 1}(d)
HilbertSchmidtStates(d::Int) = HilbertSchmidtStates{2, 1}(d)

function rand(rng::AbstractRNG, hs::HilbertSchmidtStates{β, K}) where {β, K}
    ρ = rand(rng, hs.w)
    renormalize!(ρ)
    ρ
end

#Random channels
struct ChoiJamiolkowskiMatrices{β, K} <: QIContinuousMatrixDistribution
    w::WishartEnsemble
    idim::Int
    odim::Int

    function ChoiJamiolkowskiMatrices{β, K}(idim::Int, odim::Int)  where {β, K}
        w = WishartEnsemble{β, K}(idim * odim)
        new(w, idim, odim)
    end
end

function ChoiJamiolkowskiMatrices{β}(idim::Int, odim::Int) where β
    ChoiJamiolkowskiMatrices{β, 1}(idim, odim)
end

function ChoiJamiolkowskiMatrices{β}(d::Int) where β
    ChoiJamiolkowskiMatrices{β}(d, d)
end

function ChoiJamiolkowskiMatrices(idim::Int, odim::Int)
    ChoiJamiolkowskiMatrices{2}(idim, odim)
end

function ChoiJamiolkowskiMatrices(d::Int)
    ChoiJamiolkowskiMatrices(d, d)
end

function rand(rng::AbstractRNG, c::ChoiJamiolkowskiMatrices{β, K}) where {β, K}
    z = rand(rng, c.w)
    y = ptrace(z, [c.odim, c.idim], [1])
    sy = funcmh!(x -> 1 / sqrt(x), y)
    onesy = Matrix(I, c.odim, c.odim) ⊗ sy # onesy = eye(c.odim) ⊗ sy
    DynamicalMatrix(onesy * z * onesy, c.idim, c.odim)
end

# Random POVMs implemented according to
# https://arxiv.org/pdf/1902.04751.pdf
abstract type AbstractHaarPOVM <: QIContinuousMatrixDistribution
end

struct HaarPOVM{N} <: AbstractHaarPOVM
    idim::Int
    odim::Int
    c::HaarIsometry

    function HaarPOVM{N}(idim::Int, odim::Int) where N
        c = HaarIsometry(idim::Int, N*odim::Int)
        new(idim, odim, c)
    end
end
# N controls the rank (mixedness) of the effects, N=1 gives rank-one effects
HaarPOVM(idim::Int, odim::Int) = HaarPOVM{1}(idim, odim)

#this should use slicing of V
function rand(rng::AbstractRNG, c::HaarPOVM{N}) where N
    V = rand(rng, c.c)
    POVMMeasurement([V'*(ketbra(i, i, c.odim) ⊗ 𝕀(N))*V for i=1:c.odim])
end

struct VonNeumannPOVM <: AbstractHaarPOVM
    d::Int
    c::CUE

    function VonNeumannPOVM(d::Int)
        c = CUE(d)
        new(d, c)
    end
end

function rand(rng::AbstractRNG, c::VonNeumannPOVM)
    V = rand(rng, c.c)
    POVMMeasurement([proj(V[:, i]) for i=1:c.d])
end

struct WishartPOVM{V} <: QIContinuousMatrixDistribution
    idim::Int
    odim::Int
    c::Vector{WishartEnsemble}

    function WishartPOVM{V}(idim::Int) where V
        odim = length(V)
        c = [WishartEnsemble{2, v}(idim) for v=V]
        new(idim, odim, c)
    end
end

function WishartPOVM(idim::Int, odim::Int, K::Real=1)
    V = Tuple(K .* ones(odim))
    WishartPOVM{V}(idim)
end

function rand(rng::AbstractRNG, c::WishartPOVM)
    Ws = map(x->rand(rng, x), c.c)
    S = sum(Ws)
    Ssq = funcmh!(x->1/sqrt(x), S)
    POVMMeasurement([Ssq * W * Ssq for W=Ws])
end
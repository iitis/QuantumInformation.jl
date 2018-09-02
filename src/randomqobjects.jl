using Distributions

"""
$(SIGNATURES)
- `ϕ`: vector.

Generates random ket based on `ϕ`.
"""
struct HaarKet{β} <: ContinuousMatrixDistribution
    d::Int
end

HaarKet(d::Int) = HaarKet{2}(d)

function rand(h::HaarKet{1})
    ψ = randn(h.d)
    renormalize!(ψ)
    ψ
 end

 function rand(h::HaarKet{2})
     ψ = randn(h.d) + 1im * randn(h.d)
     renormalize!(ψ)
     ψ
  end

"""
$(SIGNATURES)
- `d`: length.

Generates random ket of length `d`.
"""
struct HilbertSchmidtStates{β, K}
    w::WishartEnsemble
    d::Int

    function HilbertSchmidtStates{β, K}(d::Int) where {β, K}
        w = WishartEnsemble{β, K}(d)
        new(w, w.d)
    end
end
HilbertSchmidtStates{β}(d::Int) where β = HilbertSchmidtStates{β, 1}(d)
HilbertSchmidtStates(d::Int) = HilbertSchmidtStates{2, 1}(d)

function rand(hs::HilbertSchmidtStates{β, K}) where {β, K}
    ρ = rand(hs.w)
    renormalize!(ρ)
    ρ
end

struct ChoiJamiolkowskiMatrices{β, K}
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

function rand(c::ChoiJamiolkowskiMatrices{β, K}) where {β, K}
    z = rand(c.w)
    y = ptrace(z, [c.odim, c.idim], [1])
    sy = funcmh!(x -> 1 / sqrt(x), y)
    onesy = Matrix(I, c.odim, c.odim) ⊗ sy # onesy = eye(c.odim) ⊗ sy
    DynamicalMatrix(onesy * z * onesy, c.idim, c.odim)
end
# """
# $(SIGNATURES)
# - `ρ`: quantum state.
#
# Generates random quantum mixed state on the place of `ρ`.
# """
#
# """
# $(SIGNATURES)
# - `d`: dimension.
#
# Generates random quantum mixed state of dimension `d`.
# """
#
# """
# $(SIGNATURES)
# - `J`: dimension.
#
# Generates random Jamiolkowski state of dimension on the place of `J`.
# """
#
# """
# $(SIGNATURES)
# - `d`: dimension.
#
# Generates random Jamiolkowski state of dimension `d`.
# """

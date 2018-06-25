"""
$(SIGNATURES)
- `ϕ`: vector.

Generates random ket based on `ϕ`.
"""
function random_ket!(ϕ::AbstractVector{<:Real})
    randn!(ϕ)
    renormalize!(ϕ)
end

function random_ket!(ϕ::AbstractVector{<:Complex})
    n = length(ϕ)
    ϕ[:] = randn(n) + 1im * randn(n)
    renormalize!(ϕ)
end

"""
$(SIGNATURES)
- `d`: length.

Generates random ket of length `d`.
"""
function random_ket(::Type{T}, d::Int) where T<:Union{Real, Complex}
    ϕ=zeros(T, d)
    random_ket!(ϕ)
    ϕ
end

random_ket(d::Int) = random_ket(ComplexF64, d)

function random_mixed_state!!(ρ::AbstractMatrix{T}, A::AbstractMatrix{T}) where T<:Union{Real, Complex}
    random_ginibre_matrix!(A)
    ρ[:] = A*A'
    renormalize!(ρ)
    ρ[:] = Hermitian(ρ)
end

function random_mixed_state!(ρ::AbstractMatrix{T}, k::Int) where T<:Union{Real, Complex}
    d = size(ρ, 1)
    A = zeros(T, d, k)
    random_mixed_state!!(ρ, A)
    renormalize!(ρ)
end

"""
$(SIGNATURES)
- `ρ`: quantum state.

Generates random quantum mixed state on the place of `ρ`.
"""
function random_mixed_state!(ρ::AbstractMatrix{T}) where T<:Union{Real, Complex}
    random_mixed_state!(ρ, 1)
end

"""
$(SIGNATURES)
- `d`: dimension.

Generates random quantum mixed state of dimension `d`.
"""
function random_mixed_state(::Type{T}, d::Int64, k::Int=d) where T<:Union{Real, Complex}
    ρ = zeros(T, d, d)
    random_mixed_state!(ρ, k)
    ρ
end

random_mixed_state(d::Int64, k::Int=d) = random_mixed_state(ComplexF64, d, k)
random_mixed_state_hs(d::Int64) = random_mixed_state(d) #backwards compat

"""
$(SIGNATURES)
- `J`: dimension.

Generates random Jamiolkowski state of dimension on the place of `J`.
"""
function random_jamiolkowski_state!(J::AbstractMatrix{T}) where T<:Union{Real, Complex}
    random_dynamical_matrix!(J)
    n = isqrt(size(J, 1))
    J[:] = J[:] / n
end

"""
$(SIGNATURES)
- `d`: dimension.

Generates random Jamiolkowski state of dimension `d`.
"""
function random_jamiolkowski_state(::Type{T}, n::Int) where T<:Union{Real, Complex}
    J = zeros(T, n*n, n*n)
    random_jamiolkowski_state!(J)
    J
end

#function random_mixed_state_fixed_purity(d::Int, p::Real)
#  error("Not implemented")
#  l = random_vector_fixed_l1_l2(1., p, d)
#  u = random_unitary(d)
#  return u' * diagm(l) * u
#end

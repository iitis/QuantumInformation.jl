function random_ket!(ϕ::T) where {T<:AbstractVector{T1}} where {T1<:Real}
    randn!(ϕ)
    renormalize!(ϕ)
end

function random_ket!(ϕ::T) where {T<:AbstractVector{T1}} where {T1<:Complex}
    for i=1:length(ϕ)
        ϕ[i] = randn() + 1im * randn()
    end
    renormalize!(ϕ)
end

function random_ket(::Type{T}, d::Int) where T<:Union{Real, Complex}
    ϕ=zeros(T, d)
    random_ket!(ϕ)
    ϕ
end

random_ket(d::Int) = random_ket(ComplexF64, d)

function random_mixed_state_hs!!(ρ::AbstractMatrix{T}, A::AbstractMatrix{T}) where T<:Real
    random_ginibre_matrix!(A)
    ρ[:] = A*A'    # A_mul_Bt!(ρ, A, A) deprecated
    renormalize!(ρ)
end

function random_mixed_state_hs!!(ρ::AbstractMatrix{T}, A::AbstractMatrix{T}) where T<:Complex
    random_ginibre_matrix!(A)
    ρ[:] = A*A'    # A_mul_Bc!(ρ, A, A) deprecated
    renormalize!(ρ)
end

function random_mixed_state_hs!(ρ::AbstractMatrix{T}) where T<:Union{Real, Complex}
    A = zeros(ρ)
    random_mixed_state_hs!!(ρ, A)
    renormalize!(ρ)
end

function random_mixed_state_hs(::Type{T}, d::Int64) where T<:Union{Real, Complex}
    ρ = zeros(T, d, d)
    random_mixed_state_hs!(ρ)
    ρ
end

random_mixed_state_hs(d::Int64) = random_mixed_state_hs(ComplexF64, d)

function random_jamiolkowski_state!(J::AbstractMatrix{T}) where T<:Union{Real, Complex}
    random_dynamical_matrix!(J)
    n = round(Int, sqrt(size(J, 1)), RoundDown)
    for i=1:length(J)
        J[i] = J[i] / n
    end
end

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

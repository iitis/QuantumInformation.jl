include("randommatrix.jl")
include("utils.jl")

function random_ket!{T<:Float64}(ϕ::Vector{T})
    randn!(ϕ)
    renormalize!(ϕ)
end

function random_ket!{T<:Complex128}(ϕ::Vector{T})
    for i=1:length(ϕ)
        ϕ[i] = randn() + 1im * randn()
    end
    renormalize!(ϕ)
end

function random_ket{T<:Union(Float64, Complex128)}(M::Type{T}, d::Int64)
    ϕ=zeros(M, d)
    random_ket!(ϕ)
    ϕ
end

random_ket(d::Int64) = random_ket(Complex128, d)

function random_mixed_state_hs!!{T<:Float64}(ρ::Matrix{T}, A::Matrix{T})
    random_ginibre_matrix!(A)
    A_mul_Bt!(ρ, A, A)
    renormalize!(ρ)
end

function random_mixed_state_hs!!{T<:Complex128}(ρ::Matrix{T}, A::Matrix{T})
    random_ginibre_matrix!(A)
    A_mul_Bc!(ρ, A, A)
    renormalize!(ρ)
end

function random_mixed_state_hs!{T<:Union(Float64, Complex128)}(ρ::Matrix{T})
    A = zeros(ρ)
    random_mixed_state_hs!!(ρ, A)
    renormalize!(ρ)
end

function random_mixed_state_hs{T<:Union(Float64, Complex128)}(M::Type{T}, d::Int64)
    ρ = zeros(M, d, d)
    random_mixed_state_hs!(ρ)
    ρ
end

random_mixed_state_hs(d::Int64) = random_mixed_state_hs(Complex128, d)

function random_jamiolkowski_state!{T<:Union(Float64, Complex128)}(J::Matrix{T})
    random_dynamical_matrix!(J)
    n = int(sqrt(size(J, 1)))
    for i=1:length(J)
        J[i] = J[i] / n
    end
end

function random_jamiolkowski_state{T<:Union(Float64, Complex128)}(M::Type{T}, n::Int64)
    J = zeros(M, n*n, n*n)
    random_jamiolkowski_state!(J)
    J
end

#function random_mixed_state_fixed_purity(d::Int, p::Real)
#  error("Not implemented")
#  l = random_vector_fixed_l1_l2(1., p, d)
#  u = random_unitary(d)
#  return u' * diagm(l) * u
#end
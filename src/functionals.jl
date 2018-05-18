trace_norm(A::AbstractMatrix{T}) where T <: Number = sum(svdvals(A))

function trace_distance(A::AbstractMatrix{T}, B::AbstractMatrix{T})  where T <: Number
    one(T)/2 * trace_norm(A - B)
end

hs_norm(A::AbstractMatrix{T}) where T <: Number = sqrt(sum(svdvals(A).^2))

function hs_distance(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T <: Number
    hs_norm(A - B)
end

function fidelity_sqrt(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    throw(ArgumentError("Non square matrix"))
  end
  λ = real(eigvals(ρ * σ))
  r = sum(sqrt.(λ[λ.>0]))
end

function fidelity(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
  fidelity_sqrt(ρ, σ)^2
end

fidelity(ϕ::AbstractVector{T}, ψ::AbstractVector{T}) where T<:Number = abs2(dot(ϕ, ψ))
fidelity(ϕ::AbstractVector{T}, ρ::AbstractMatrix{T}) where T<:Number = ϕ' * ρ * ϕ
fidelity(ρ::AbstractMatrix{T}, ϕ::AbstractVector{T}) where T<:Number = fidelity(ϕ, ρ)

function gate_fidelity(U::AbstractMatrix{T}, V::AbstractMatrix{T}) where T<:Number
    abs(1.0 / size(U,1) * trace(U'*V))
end

shannon_entropy(p::AbstractVector{T}) where T<:Real = -sum(p .* log.(p))
shannon_entropy(x::T) where T<:Real = x > 0 ? -x * log(x) - (1 - x) * log(1 - x) : error("Negative number passed to shannon_entropy")

"""
Calculates the von Neuman entropy of positive matrix \rho
S(\rho)=-tr(\rho\log(\rho))

Equivalent faster form:
S(\rho)=-\sum_i \lambda_i(\rho)*\log(\lambda_i(\rho))
http://www.quantiki.org/wiki/Von_Neumann_entropy
"""
function quantum_entropy(ρ::Hermitian{T}) where T<:Number
    λ = eigvals(ρ)
    shannon_entropy(λ[λ .> 0])
end

quantum_entropy(H::AbstractMatrix{T}) where T<:Number = ishermitian(H) ? quantum_entropy(Hermitian(H)) : error("Non-hermitian matrix passed to entropy")
quantum_entropy(ϕ::AbstractVector{T}) where T<:Number = zero(T)

function relative_entropy(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    log_σ = funcmh(log, σ)
    real(-quantum_entropy(ρ) - trace(ρ * log_σ))
end

function kl_divergence(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    relative_entropy(ρ, σ)
end

function js_divergence(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    0.5kl_divergence(ρ, σ) + 0.5kl_divergence(σ, ρ)
end

"""
D_B(\rho,\sigma) = \sqrt{2-2\sqrt{F(\rho,\sigma)}}

http://www.quantiki.org/wiki/Fidelity
"""
function bures_distance(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    sqrt(2 - 2 * fidelity_sqrt(ρ, σ))
end

"""
D_A(\rho,\sigma) = \arccos \sqrt{F(\rho,\sigma)}

http://www.quantiki.org/wiki/Fidelity
"""
function bures_angle(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    acos(fidelity_sqrt(ρ, σ))
end

"""
G(\rho,\sigma) = \mathrm{tr}\rho\sigma + \sqrt{1-\mathrm{tr}(\rho^2)}\sqrt{1-\mathrm{tr}(\sigma^2)}

http://www.quantiki.org/wiki/Superfidelity
``Sub-- and super--fidelity as bounds for quantum fidelity''
Quantum Information & Computation, Vol.9 No.1&2 (2009)
"""
function superfidelity(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    return trace(ρ'*σ) + sqrt(1 - trace(ρ'*ρ)) * sqrt(1 - trace(σ'*σ))
end

"""
    Implement according to
    title={Quantum entanglement},
author={Horodecki, R. and Horodecki, P. and Horodecki, M. and Horodecki, K.},
journal={Reviews of Modern Physics},
"""
function negativity(ρ::AbstractMatrix, dims::Vector, sys::Int)
    ρ_s = ptranspose(ρ, dims, sys)
    λ = eigvals(ρ_s)
    -real(sum(λ[λ .< 0]))
end

function log_negativity(ρ::AbstractMatrix, dims::Vector, sys::Int)
    ρ_s = ptranspose(ρ, dims, sys)
    log(trace_norm(ρ_s))
end

function ppt(ρ::AbstractMatrix, dims::Vector, sys::Int)
    ρ_s = ptranspose(ρ, dims, sys)
    minimum(eigvals(ρ_s))
end

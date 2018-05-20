using Convex, SCS

trace_norm(A::AbstractMatrix{T}) where T <: Number = sum(svdvals(A))

function trace_distance(A::AbstractMatrix{T}, B::AbstractMatrix{T})  where T <: Number
    one(T)/2 * trace_norm(A - B)
end

hs_norm(A::AbstractMatrix{T}) where T <: Number = sqrt(sum(abs2.(A)))

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

function bures_distance(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    sqrt(2 - 2 * fidelity_sqrt(ρ, σ))
end

function bures_angle(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    acos(fidelity_sqrt(ρ, σ))
end


function superfidelity(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    return trace(ρ'*σ) + sqrt(1 - trace(ρ'*ρ)) * sqrt(1 - trace(σ'*σ))
end

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

function diamond_norm(J::AbstractMatrix{T}, d1::Int, d2::Int=d1) where T <: Number
  X = ComplexVariable(d1*d2, d1*d2)
  t = 0.5inner_product(X, J)+0.5inner_product(X', J')

  ρ₀ = ComplexVariable(d1, d1)
  ρ₁ = ComplexVariable(d1, d1)

  constraints = [ρ₀ in :SDP, ρ₁ in :SDP]
  constraints += trace(ρ₀) == 1
  constraints += trace(ρ₁) == 1
  constraints += [eye(d2) ⊗ ρ₀ X; X' eye(d2) ⊗ ρ₁] in :SDP

  problem = maximize(t, constraints)
  solve!(problem, SCSSolver(verbose=0))
  problem.optval
end

function diamond_distance(J1::AbstractMatrix{T}, J2::AbstractMatrix{T}, d1::Int, d2::Int=d1) where T <: Number
    diamond_norm(J1-J2, d1, d2)
end

using Convex, SCS

"""
$(SIGNATURES)
- `A`: matrix.

Return [trace norm](https://www.quantiki.org/wiki/trace-norm) of matrix `A`.
"""
trace_norm(A::AbstractMatrix{T}) where T <: Number = sum(svdvals(A))

"""
$(SIGNATURES)
- `A`: matrix.
- `B`: matrix.

Return [trace distance](https://www.quantiki.org/wiki/trace-distance) between matrices `A` and `B`.
"""
function trace_distance(A::AbstractMatrix{T}, B::AbstractMatrix{T})  where T <: Number
    one(T)/2 * trace_norm(A - B)
end

"""
$(SIGNATURES)
- `A`: matrix.

Return [Hilbert–Schmidt norm](https://en.wikipedia.org/wiki/Hilbert%E2%80%93Schmidt_operator) of matrix `A`.
"""
hs_norm(A::AbstractMatrix{T}) where T <: Number = sqrt(sum(abs2.(A)))

"""
$(SIGNATURES)
- `A`: matrix.
- `B`: matrix.

Return [Hilbert–Schmidt distance](https://en.wikipedia.org/wiki/Hilbert%E2%80%93Schmidt_operator) between matrices `A` and `B`.
"""
function hs_distance(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T <: Number
    hs_norm(A - B)
end

"""
$(SIGNATURES)
- `ρ`: matrix.
- `σ`: matrix.

Return square root of [fidelity](https://www.quantiki.org/wiki/fidelity) between matrices `ρ` and `σ`.
"""
function fidelity_sqrt(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    throw(ArgumentError("Non square matrix"))
  end
  λ = real(eigvals(ρ * σ))
  r = sum(sqrt.(λ[λ.>0]))
end

"""
$(SIGNATURES)
- `ρ`: matrix.
- `σ`: matrix.

Return [fidelity](https://www.quantiki.org/wiki/fidelity) between matrices `ρ` and `σ`.
"""
function fidelity(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
  fidelity_sqrt(ρ, σ)^2
end

fidelity(ϕ::AbstractVector{T}, ψ::AbstractVector{T}) where T<:Number = abs2(dot(ϕ, ψ))
fidelity(ϕ::AbstractVector{T}, ρ::AbstractMatrix{T}) where T<:Number = ϕ' * ρ * ϕ
fidelity(ρ::AbstractMatrix{T}, ϕ::AbstractVector{T}) where T<:Number = fidelity(ϕ, ρ)

"""
$(SIGNATURES)
- `U`: quantum gate.
- `V`: quantum gate.

Return [fidelity](https://www.quantiki.org/wiki/fidelity) between gates `U` and `V`.
"""
function gate_fidelity(U::AbstractMatrix{T}, V::AbstractMatrix{T}) where T<:Number
    abs(1.0 / size(U,1) * trace(U'*V))
end

"""
$(SIGNATURES)
- `p`: vector.

Return [Shannon entorpy](https://en.wikipedia.org/wiki/Entropy_(information_theory)) of vector `p`.
"""
shannon_entropy(p::AbstractVector{T}) where T<:Real = -sum(p .* log.(p))

"""
$(SIGNATURES)
- `x`: real number.

Return binary [Shannon entorpy](https://en.wikipedia.org/wiki/Entropy_(information_theory)) given by \$-x  \\log(x) - (1 - x)  \\log(1 - x)\$.
"""
shannon_entropy(x::T) where T<:Real = x > 0 ? -x * log(x) - (1 - x) * log(1 - x) : error("Negative number passed to shannon_entropy")

"""
$(SIGNATURES)
- `ρ`: quantum state.

Return [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy) of quantum state `ρ`.
"""
function quantum_entropy(ρ::Hermitian{T}) where T<:Number
    λ = eigvals(ρ)
    shannon_entropy(λ[λ .> 0])
end

quantum_entropy(H::AbstractMatrix{T}) where T<:Number = ishermitian(H) ? quantum_entropy(Hermitian(H)) : error("Non-hermitian matrix passed to entropy")
quantum_entropy(ϕ::AbstractVector{T}) where T<:Number = zero(T)

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [quantum relative entropy](https://en.wikipedia.org/wiki/Quantum_relative_entropy) of quantum state `ρ` with respect to `σ`.
"""
function relative_entropy(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    log_σ = funcmh(log, σ)
    real(-quantum_entropy(ρ) - trace(ρ * log_σ))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [Kullback–Leibler divergence](https://en.wikipedia.org/wiki/Quantum_relative_entropy) of quantum state `ρ` with respect to `σ`.
"""
function kl_divergence(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    relative_entropy(ρ, σ)
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [Jensen–Shannon divergence](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence) of quantum state `ρ` with respect to `σ`.
"""
function js_divergence(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    0.5kl_divergence(ρ, σ) + 0.5kl_divergence(σ, ρ)
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [Bures distance](https://en.wikipedia.org/wiki/Bures_metric) between quantum states `ρ` and `σ`.
"""
function bures_distance(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    sqrt(2 - 2 * fidelity_sqrt(ρ, σ))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [Bures angle](https://en.wikipedia.org/wiki/Bures_metric) between quantum states `ρ` and `σ`.
"""
function bures_angle(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    acos(fidelity_sqrt(ρ, σ))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [superfidelity](https://www.quantiki.org/wiki/superfidelity) between quantum states `ρ` and `σ`.
"""
function superfidelity(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    return trace(ρ'*σ) + sqrt(1 - trace(ρ'*ρ)) * sqrt(1 - trace(σ'*σ))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `dims`: dimensions of subsystems.
- `sys`: transposed subsystem.

Return [negativity](https://www.quantiki.org/wiki/negativity) of quantum state `ρ`.
"""
function negativity(ρ::AbstractMatrix, dims::Vector, sys::Int)
    ρ_s = ptranspose(ρ, dims, sys)
    λ = eigvals(ρ_s)
    -real(sum(λ[λ .< 0]))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `dims`: dimensions of subsystems.
- `sys`: transposed subsystem.

Return [log negativity](https://www.quantiki.org/wiki/negativity) of quantum state `ρ`.
"""
function log_negativity(ρ::AbstractMatrix, dims::Vector, sys::Int)
    ρ_s = ptranspose(ρ, dims, sys)
    log(trace_norm(ρ_s))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `dims`: dimensions of subsystems.
- `sys`: transposed subsystem.

Return minimum eigenvalue of [positive partial transposition](https://www.quantiki.org/wiki/positive-partial-transpose) of quantum state `ρ`.
"""
function ppt(ρ::AbstractMatrix, dims::Vector, sys::Int)
    ρ_s = ptranspose(ρ, dims, sys)
    minimum(eigvals(ρ_s))
end

"""
$(SIGNATURES)
- `J`: matrix.
- `d1`: ?.
- `d2`: ?.

Return [diamond norm](https://arxiv.org/pdf/1004.4110.pdf) of matrix `J`.
"""
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

"""
$(SIGNATURES)
- `J1`: matrix.
- `J2`: matrix.
- `d1`: ?.
- `d2`: ?.

Return [diamond distance](https://arxiv.org/pdf/1004.4110.pdf) between matrices `J1` and `J2`.
"""
function diamond_distance(J1::AbstractMatrix{T}, J2::AbstractMatrix{T}, d1::Int, d2::Int=d1) where T <: Number
    diamond_norm(J1-J2, d1, d2)
end

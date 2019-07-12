export norm_trace, trace_distance, norm_hs, hs_distance,
    purity, fidelity_sqrt, fidelity, gate_fidelity, shannon_entropy,
    vonneumann_entropy, renyi_entropy, relative_entropy, kl_divergence,
    js_divergence, bures_distance, bures_angle, superfidelity, negativity,
    log_negativity, ppt, concurrence

"""
$(SIGNATURES)
- `A`: matrix.

Return [trace norm](https://www.quantiki.org/wiki/trace-norm) of matrix `A`.
"""
norm_trace(A::AbstractMatrix{<:Number}) = sum(svdvals(A))

"""
$(SIGNATURES)
- `A`: matrix.
- `B`: matrix.

Return [trace distance](https://www.quantiki.org/wiki/trace-distance) between matrices `A` and `B`.
"""
function trace_distance(A::AbstractMatrix{T1}, B::AbstractMatrix{T2}) where {T1<:Number, T2<:Number}
    norm_trace(A - B) / 2
end

"""
$(SIGNATURES)
- `A`: matrix.

Return [Hilbert–Schmidt norm](https://en.wikipedia.org/wiki/Hilbert%E2%80%93Schmidt_operator) of matrix `A`.
"""
norm_hs(A::AbstractMatrix{<:Number}) = sqrt(sum(abs2.(A)))

"""
$(SIGNATURES)
- `A`: matrix.
- `B`: matrix.

Return [Hilbert–Schmidt distance](https://en.wikipedia.org/wiki/Hilbert%E2%80%93Schmidt_operator) between matrices `A` and `B`.
"""
function hs_distance(A::AbstractMatrix{<:Number}, B::AbstractMatrix{<:Number})
    norm_hs(A - B)
end

"""
$(SIGNATURES)
- `ρ`: matrix.

Return the purity of `ρ` ∈ [1/d, 1]
"""
purity(ρ::AbstractMatrix{<:Number}) = tr(ρ^2)

"""
$(SIGNATURES)
- `ρ`: matrix.
- `σ`: matrix.

Return square root of [fidelity](https://www.quantiki.org/wiki/fidelity) between matrices `ρ` and `σ`.
"""
function fidelity_sqrt(ρ::AbstractMatrix{<:Number}, σ::AbstractMatrix{<:Number})
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    throw(ArgumentError("Non square matrix"))
  end
  λ = real(eigvals(ρ * σ))
  r = real(sum(sqrt.(λ[λ.>0])))
end

"""
$(SIGNATURES)
- `ρ`: matrix.
- `σ`: matrix.

Return [fidelity](https://www.quantiki.org/wiki/fidelity) between matrices `ρ` and `σ`.
"""
function fidelity(ρ::AbstractMatrix{<:Number}, σ::AbstractMatrix{<:Number})
  fidelity_sqrt(ρ, σ)^2
end

fidelity(ϕ::AbstractVector{<:Number}, ψ::AbstractVector{<:Number}) = abs2(dot(ϕ, ψ))
fidelity(ϕ::AbstractVector{<:Number}, ρ::AbstractMatrix{<:Number}) = real(ϕ' * ρ * ϕ)
fidelity(ρ::AbstractMatrix{<:Number}, ϕ::AbstractVector{<:Number}) = fidelity(ϕ, ρ)

"""
$(SIGNATURES)
- `U`: quantum gate.
- `V`: quantum gate.

Return [fidelity](https://www.quantiki.org/wiki/fidelity) between gates `U` and `V`.
"""
function gate_fidelity(U::AbstractMatrix{<:Number}, V::AbstractMatrix{<:Number})
    abs(1.0 / size(U,1) * tr(U'*V))
end

"""
$(SIGNATURES)
- `p`: vector.

Return [Shannon entorpy](https://en.wikipedia.org/wiki/Entropy_(information_theory)) of vector `p`.
"""
shannon_entropy(p::AbstractVector{<:Real}) = -sum(p .* log.(p))

"""
$(SIGNATURES)
- `x`: real number.

Return binary [Shannon entorpy](https://en.wikipedia.org/wiki/Entropy_(information_theory)) given by \$-x  \\log(x) - (1 - x)  \\log(1 - x)\$.
"""
shannon_entropy(x::Real) = x > 0 ? -x * log(x) - (1 - x) * log(1 - x) : error("Negative number passed to shannon_entropy")

"""
$(SIGNATURES)
- `ρ`: quantum state.

Return [Von Neumann entropy](https://en.wikipedia.org/wiki/Von_Neumann_entropy) of quantum state `ρ`.
"""
function vonneumann_entropy(ρ::Hermitian{<:Number})
    λ = eigvals(ρ)
    shannon_entropy(λ[λ .> 0])
end

vonneumann_entropy(H::AbstractMatrix{<:Number}) = ishermitian(H) ? vonneumann_entropy(Hermitian(H)) : error("Non-hermitian matrix passed to entropy")
vonneumann_entropy(ϕ::AbstractVector{T}) where T<:Number = zero(T)

function renyi_entropy(ρ::Hermitian{<:Number}, α::Real)
    α >= 0 && α!=1 ? () : throw(ArgumentError("Parameter α must be α ≥ 0 and α ≠ 1"))
    λ = eigvals(ρ)
    λ = λ[λ .> 0]
    1/(1-α)*log(sum(λ.^α))
end

renyi_entropy(ρ::AbstractMatrix{<:Number}, α::Real) = renyi_entropy(Hermitian(ρ), α)
"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [quantum relative entropy](https://en.wikipedia.org/wiki/Quantum_relative_entropy) of quantum state `ρ` with respect to `σ`.
"""
function relative_entropy(ρ::AbstractMatrix{<:Number}, σ::AbstractMatrix{<:Number})
    log_σ = funcmh(log, σ)
    real(-vonneumann_entropy(ρ) - tr(ρ * log_σ))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [Kullback–Leibler divergence](https://en.wikipedia.org/wiki/Quantum_relative_entropy) of quantum state `ρ` with respect to `σ`.
"""
function kl_divergence(ρ::AbstractMatrix{<:Number}, σ::AbstractMatrix{<:Number})
    relative_entropy(ρ, σ)
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [Jensen–Shannon divergence](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence) of quantum state `ρ` with respect to `σ`.
"""
function js_divergence(ρ::AbstractMatrix{<:Number}, σ::AbstractMatrix{<:Number})
    0.5kl_divergence(ρ, σ) + 0.5kl_divergence(σ, ρ)
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [Bures distance](https://en.wikipedia.org/wiki/Bures_metric) between quantum states `ρ` and `σ`.
"""
function bures_distance(ρ::AbstractMatrix{<:Number}, σ::AbstractMatrix{<:Number})
    sqrt(2 - 2 * fidelity_sqrt(ρ, σ))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [Bures angle](https://en.wikipedia.org/wiki/Bures_metric) between quantum states `ρ` and `σ`.
"""
function bures_angle(ρ::AbstractMatrix{<:Number}, σ::AbstractMatrix{<:Number})
    acos(fidelity_sqrt(ρ, σ))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `σ`: quantum state.

Return [superfidelity](https://www.quantiki.org/wiki/superfidelity) between quantum states `ρ` and `σ`.
"""
function superfidelity(ρ::AbstractMatrix{<:Number}, σ::AbstractMatrix{<:Number})
    return real(tr(ρ'*σ) + sqrt(1 - tr(ρ'*ρ)) * sqrt(1 - tr(σ'*σ)))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `dims`: dimensions of subsystems.
- `sys`: transposed subsystem.

Return [negativity](https://www.quantiki.org/wiki/negativity) of quantum state `ρ`.
"""
function negativity(ρ::AbstractMatrix{<:Number}, dims::Vector{Int}, sys::Int)
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
function log_negativity(ρ::AbstractMatrix{<:Number}, dims::Vector{Int}, sys::Int)
    ρ_s = ptranspose(ρ, dims, sys)
    log(norm_trace(ρ_s))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `dims`: dimensions of subsystems.
- `sys`: transposed subsystem.

Return minimum eigenvalue of [positive partial transposition](https://www.quantiki.org/wiki/positive-partial-transpose) of quantum state `ρ`.
"""
function ppt(ρ::AbstractMatrix{<:Number}, dims::Vector{Int}, sys::Int)
    ρ_s = ptranspose(ρ, dims, sys)
    minimum(eigvals(ρ_s))
end

"""
$(SIGNATURES)
- `ρ`: quantum state.

Calculates the [concurrence of a two-qubit system](https://en.wikipedia.org/wiki/Concurrence_(quantum_computing)) `ρ`.
"""
function concurrence(ρ::AbstractMatrix{<:Number})
    if size(ρ, 1) != 4
        throw(ArgumentError("Concurrence only properly defined for two qubit systems"))
    end

    σ = (sy ⊗ sy)*conj.(ρ)*(sy ⊗ sy)
    λ = sqrt.(sort(real(eigvals(ρ*σ)), rev=true))
    max(0, λ[1]-λ[2]-λ[3]-λ[4])
end

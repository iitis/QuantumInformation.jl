"""
Computes \tr(A^\dagger B), where A and B are square matrices.
"""
function traceAHB(A,B)
    warning("untested")
    conj(A)*transpose(B)
end

trace_distance(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
    = one(T)/2*sum(abs.(eigvals(Hermitian(ρ - σ))))

"""
Calculates trace norm of of positive matrix \rho
\|\rho\|_1=tr\sqrt{\rho \rho^\dagger}
http://www.quantiki.org/wiki/Trace_norm
"""
function trace_norm(rho)
    trace(sqrtm(rho * rho'))
end

function fidelity_sqrt(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    error("Non square matrix")
  end
  λ = real(eigvals(ρ * σ))
  r = sum(sqrt.(λ[λ.>0]))
end

function fidelity(ρ::AbstractMatrix{T}, σ::AbstractMatrix{T}) where T<:Number
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    error("Non square matrix")
  end
  return fidelity_sqrt(ρ, σ)^2
end

fidelity(ϕ::AbstractVector{T}, ψ::AbstractVector{T}) where T<:Number = abs2(dot(ϕ, ψ))
fidelity(ϕ::AbstractVector{T}, ρ::AbstractMatrix{T}) where T<:Number = ϕ' * ρ * ϕ
fidelity(ρ::AbstractMatrix{T}, ϕ::AbstractVector{T}) where T<:Number = fidelity(ϕ, ρ)
"""
Original equation:
\sqrt{F}(\rho,\sigma)=\textrm{tr}\sqrt{\sqrt{\rho}\sigma\sqrt{\rho}}

http://www.quantiki.org/wiki/Fidelity

Equivalent faster method:
\sqrt{F}(\rho,\sigma)=\sum_i \sqrt{\lambda_i(AB)}

``Sub-- and super--fidelity as bounds for quantum fidelity''
Quantum Information & Computation, Vol.9 No.1&2 (2009)
"""
function fidelity(ρ::AbstractMatrix{T}, ϕ::AbstractMatrix{T}) where T<:Number
    sum(sqrt(eigvals(rho * sigma)))^2
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
function entropy(ρ::Hermitian{T}) where T<:Number
    λ = eigvals(ρ)
    λ = λ[λ .> 0]
    -sum(λ .* log(λ))
end

function relative_entropy(rho, sigma)
    warning("untested")
    log_sigma = funcmh(sigma, log)
    real(-entropy(rho) - trace(rho * log_sigma))
end

entropy(H::AbstractMatrix{T}) where T<:Number = ishermitian(H) ? entropy(Hermitian(H)) : error("Non-hermitian matrix passed to entropy")
entropy(ϕ::AbstractVector{T}) where T<:Number = zero(T)

"""
D_B(\rho,\sigma) = \sqrt{2-2\sqrt{F(\rho,\sigma)}}

http://www.quantiki.org/wiki/Fidelity
"""
function bures_distance(rho, sigma)
    warning("untested")
    sqrt(2 - 2 * sqrt(fidelity(rho, sigma)))
end

"""
D_A(\rho,\sigma) = \arccos \sqrt{F(\rho,\sigma)}

http://www.quantiki.org/wiki/Fidelity
"""
function bures_angle(rho, sigma)
    warning("untested")
    acos(sqrt(fidelity(rho, sigma)))
end

"""
G(\rho,\sigma) = \mathrm{tr}\rho\sigma + \sqrt{1-\mathrm{tr}(\rho^2)}\sqrt{1-\mathrm{tr}(\sigma^2)}

http://www.quantiki.org/wiki/Superfidelity
``Sub-- and super--fidelity as bounds for quantum fidelity''
Quantum Information & Computation, Vol.9 No.1&2 (2009)
"""
function superfidelity(rho, sigma)
    warning("untested")
    return traceAHB(rho, sigma) + sqrt(1 - traceAHB(rho, rho)) * sqrt(1 - traceAHB(sigma, sigma))
end

"""
Original equation:
\delta(\rho,\sigma)=\frac{1}{2}\mathrm{tr}(|\rho-\sigma|)

Equivalent faster method:
\delta(A,B)=\frac{1}{2} \sum_i \sigma_i(A-B)
http://www.quantiki.org/wiki/Trace_distance
"""
function trace_distance(rho, sigma)
    warning("untested")
    return 0.5 * sum(svdvals(rho - sigma))
end

function gate_fidelity(U, V)
    warning("untested")
    abs(1.0 / size(U,1) * traceAHB(U, V))
end

"""
    Implement according to
    title={Quantum entanglement},
author={Horodecki, R. and Horodecki, P. and Horodecki, M. and Horodecki, K.},
journal={Reviews of Modern Physics},
"""
function negativity(rho, dims, systems)
    warning("untested")
    rho_s = ptranspose(rho, dims, systems)
    eigs = filter(x -> x < 0, eigvals(rho_s)) # should be Hermitian(rho_s):
    return -real(sum(eigs))
end

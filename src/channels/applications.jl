export applychannel

################################################################################
# Application of channels
################################################################################
"""
  - `Φ`: dynamical matrix.
  - `ρ`: quantum state.

Application of dynamical matrix into state `ρ`.
"""
function applychannel(
    Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}},
    ρ::AbstractMatrix{T},
) where {T <: Number}
    return ptrace(Φ.matrix * (Diagonal{T}(I, Φ.odim)⊗transpose(ρ)), [Φ.odim, Φ.idim], [2])
end

"""
  - `Φ`: list of vectors.
  - `ρ`: input matrix.

Return application of channel `Φ`` on `ρ`. Kraus representation of quantum channel
\$\\Phi\$ is a set \$\\{K_i\\}_{i\\in I}\$ of bounded operators on \$\\mathcal{H}\$
such that \$\\sum_{i\\in I} K_i^\\dagger K_i = \\mathcal{1}\$.
Then \$\\Phi(\\rho)=\\sum_{i\\in I} K_i \\rho K_i^\\dagger\$.
"""
function applychannel(
    Φ::KrausOperators{<:AbstractMatrix{<:Number}},
    ρ::AbstractMatrix{<:Number},
)
    return sum(k * ρ * k' for k in Φ.matrices)
end

"""
  - `Φ`: super-operator matrix.
  - `ρ`: quantum state.

Application of super-operator matrix into state `ρ`.
"""
function applychannel(
    Φ::SuperOperator{<:AbstractMatrix{<:Number}},
    ρ::AbstractMatrix{<:Number},
)
    return unres(Φ.matrix * res(ρ))
end

"""
  - `Φ`: Stinespring representation of quantum channel.
  - `ρ`: quantum state.
  - `dims`: dimensions of registers of `ρ`.

Application of Stinespring representation of quantum channel into state `ρ`.
"""
function applychannel(
    Φ::Stinespring{<:AbstractMatrix{<:Number}},
    ρ::AbstractMatrix{<:Number},
)
    s = Φ.matrix * ρ * Φ.matrix'
    return ptrace(s, [Φ.odim, Φ.odim*Φ.idim], [2])
end

"""
  - `Φ`: Identity channel.
  - `ρ`: quantum state.

Return application of Identity channel `Φ` on `ρ`.
"""
function applychannel(
    Φ::IdentityChannel{<:AbstractMatrix{<:Number}},
    ρ::AbstractMatrix{<:Number},
)
    # TODO: promote type
    return ρ
end

"""
  - `Φ`: Unitary channel.
  - `ρ`: quantum state.

Return application of Unitary channel `Φ` on `ρ`.
"""
function applychannel(
    Φ::UnitaryChannel{<:AbstractMatrix{<:Number}},
    ρ::AbstractMatrix{<:Number},
)
    # TODO: promote type
    return Φ.matrix*ρ*Φ.matrix'
end

"""
  - `Φ`: quantum channel.
  - `ψ`: quantum state vector.

Return application of channel `Φ` on state vector `ψ`.
"""
function applychannel(Φ::AbstractQuantumOperation, ψ::AbstractVector{<:Number})
    # TODO: promote type
    return applychannel(Φ, proj(ψ))
end

"""
  - `Φ`: Unitary channel.
  - `ψ`: quantum state vector.

Return application of Unitary channel `Φ` on state vector `ψ`.
"""
function applychannel(
    Φ::UnitaryChannel{<:AbstractMatrix{<:Number}},
    ψ::AbstractVector{<:Number},
)
    # TODO: promote type
    return Φ.matrix*ψ
end

"""
  - `Φ`: Identity channel.
  - `ψ`: quantum state vector.

Return application of Identity channel `Φ` on state vector `ψ`.
"""
function applychannel(
    Φ::IdentityChannel{<:AbstractMatrix{<:Number}},
    ψ::AbstractVector{<:Number},
)
    # TODO: promote type
    return ψ
end

function applychannel(
    Φ::POVMMeasurement{T},
    ρ::AbstractMatrix{<:Number},
) where {T <: AbstractMatrix{<:Number}}
    # TODO: Check if idim and odim are compatible with matrix and length
    probs = [real(tr(p'*ρ)) for p in Φ.matrices]
    return Diagonal(probs)
end

function applychannel(
    Φ::PostSelectionMeasurement{T},
    ρ::AbstractMatrix{<:Number},
) where {T <: AbstractMatrix{<:Number}}
    # TODO: Check if idim and odim are compatible with matrix and length
    e = Φ.matrix
    return e * ρ * e'
end

################################################################################
# making channels callable
################################################################################
for qop in (
    :KrausOperators,
    :SuperOperator,
    :DynamicalMatrix,
    :Stinespring,
    :UnitaryChannel,
    :IdentityChannel,
    :POVMMeasurement,
    :PostSelectionMeasurement,
)
    @eval begin
        function (Φ::$qop)(ρ)
            return applychannel(Φ, ρ)
        end
    end
end

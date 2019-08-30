export applychannel

################################################################################
# Application of channels
################################################################################
"""
$(SIGNATURES)
- `Φ`: dynamical matrix.
- `ρ`: quantum state.

Application of dynamical matrix into state `ρ`.
"""
function applychannel(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{T}) where T<:Number
    ptrace(Φ.matrix * (Diagonal{T}(I, Φ.odim)⊗transpose(ρ)), [Φ.odim, Φ.idim], [2])
end

"""
$(SIGNATURES)
- `Φ`: list of vectors.
- `ρ`: input matrix.

Return application of channel `Φ`` on `ρ`. Kraus representation of quantum channel
\$\\Phi\$ is a set \$\\{K_i\\}_{i\\in I}\$ of bounded operators on \$\\mathcal{H}\$
such that \$\\sum_{i\\in I} K_i^\\dagger K_i = \\mathcal{1}\$.
Then \$\\Phi(\\rho)=\\sum_{i\\in I} K_i \\rho K_i^\\dagger\$.
"""
function applychannel(Φ::KrausOperators{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    sum(k * ρ * k' for k in Φ.matrices)
end

"""
$(SIGNATURES)
- `Φ`: super-operator matrix.
- `ρ`: quantum state.

Application of super-operator matrix into state `ρ`.
"""
function applychannel(Φ::SuperOperator{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    unres(Φ.matrix * res(ρ))
end

"""
$(SIGNATURES)
- `Φ`: Stinespring representation of quantum channel.
- `ρ`: quantum state.
- `dims`: dimensions of registers of `ρ`.

Application of Stinespring representation of quantum channel into state `ρ`.
"""
function applychannel(Φ::Stinespring{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    s = Φ.matrix * ρ * Φ.matrix'
    ptrace(s, [Φ.odim, Φ.odim*Φ.idim], [2])
end

function applychannel(Φ::IdentityChannel{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    # TODO: promote type
    ρ
end

function applychannel(Φ::UnitaryChannel{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    # TODO: promote type
    Φ.matrix*ρ*Φ.matrix'
end

function applychannel(Φ::AbstractQuantumOperation, ψ::AbstractVector{<:Number})
    # TODO: promote type
    applychannel(Φ, proj(ψ))
end

function applychannel(Φ::UnitaryChannel{<:AbstractMatrix{<:Number}}, ψ::AbstractVector{<:Number})
    # TODO: promote type
    Φ.matrix*ψ
end

function applychannel(Φ::IdentityChannel{<:AbstractMatrix{<:Number}}, ψ::AbstractVector{<:Number})
    # TODO: promote type
    ψ
end

function applychannel(Φ::POVMMeasurement{T}, ρ::AbstractMatrix{<:Number}) where T<:AbstractMatrix{<:Number}
    # TODO: Check if idim and odim are compatible with matrix and length
    probs = [real(tr(p'*ρ)) for p in Φ.matrices]
    Diagonal(probs)
end

function applychannel(Φ::PostSelectionMeasurement{T}, ρ::AbstractMatrix{<:Number}) where T<:AbstractMatrix{<:Number}
    # TODO: Check if idim and odim are compatible with matrix and length
    e = Φ.matrix
    e * ρ * e'
end

################################################################################
# making channels callable
################################################################################
for qop in (:KrausOperators, :SuperOperator, :DynamicalMatrix, :Stinespring,
    :UnitaryChannel, :POVMMeasurement, :PostSelectionMeasurement)
    @eval begin
        function (Φ::$qop)(ρ)
            applychannel(Φ, ρ)
        end
    end
end

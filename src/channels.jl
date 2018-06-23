import Base.convert
import Base.size
import Base.length
import Base.kron

abstract type AbstractQuantumOperation{T} end

struct KrausOperators{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrices::Vector{T}
    idim::Int
    odim::Int
    function KrausOperators(v::Vector{T}) where T<:AbstractMatrix{<:Number}
        sizes = [size(k) for k in v]
        for s in sizes[2:end]
            if s!=sizes[1]
                throw(ArgumentError("Kraus operators list contains matrices of different dimmension"))
            end
        end
        idim, odim = sizes[1]
        new{T}(v, idim, odim)
    end
end

length(Φ::KrausOperators) = length(Φ.matrices)
size(Φ::KrausOperators) = (idim, odim)

# TODO: create iterator over KrausOperators see: https://docs.julialang.org/en/v0.6.3/manual/interfaces/

struct SuperOperator{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrix::T
    idim::Int
    odim::Int
    function SuperOperator(m::T) where T<:AbstractMatrix{<:Number}
        r, c = sizes(m)
        sr = isqrt(r), sc = isqrt(c)
        if r!=sr^2 || c!=sc^2
            throw(ArgumentError("Superoperator matrix has bad dimensions"))
        end
        idim, odim = sr, sc
        new{T}(v, idim, odim)
    end
end

struct DynamicalMatrix{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrix::T
    idim::Int
    odim::Int
    # TODO: write inner constructor
end

struct Stinespring{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrix::T
    idim::Int
    odim::Int
    # TODO: write inner constructor
end

struct IdentityChannel{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    idim::Int
    odim::Int
end

for qop in (:SuperOperator, :DynamicalMatrix, :Stinespring)
    @eval size(Φ::$qop) = size(Φ.matrix)
end

for qop in (:KrausOperators, :SuperOperator, :DynamicalMatrix, :Stinespring)
    @eval begin
        function (Φ::$qop)(ρ)
            applychannel(Φ, ρ)
        end
    end
end

for qop in (:SuperOperator, :DynamicalMatrix, :Stinespring)
    @eval convert(::Type{<:AbstractMatrix{<:Number}}, Φ::$qop) = Φ.matrix
end


"""_
$(SIGNATURES)
- `kraus`: list of Kraus operators.
- `atol`: tolerance of approximation.

Checks if set of Kraus operators fulfill completness relation.
"""
function iscptp(Φ::KrausOperators{<:AbstractMatrix{<:Number}}, atol=1e-08)
    complentess_relation = sum(k'*k for k in Φ.matrices)
    isapprox(complentess_relation, eye(complentess_relation), atol=atol)
end

function iscptp(s::SuperOperator)
    #TODO: implement
end

function iscptp(s::DynamicalMatrix)
    #TODO: implement
end

function iscptp(s::Stinespring)
    #TODO: implement
end

#### Relationship among representations of channels

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into super-operator matrix.
"""
function Base.convert(::Type{SuperOperator{T}}, Φ::KrausOperators{T}) where T<:AbstractMatrix{<:Number}
    SuperOperator(sum(k⊗(conj.(k)) for k in Φ.matrices), Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- `channel`: quantum channel map.
- `dim`: square root of the super-operator matrix dimension.

Transforms quntum channel into super-operator matrix.
"""
function SuperOperator{T}(channel::Function, idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
    odim > 0 ? () : error("Channel dimension has to be nonnegative") # TODO: fix

    m = zeros(T, idim^2, odim^2)
    for (i, e) in enumerate(base_matrices(idim)) # TODO : base_matrices should be not only square
        m[:, i] = res(channel(e))
    end
    SuperOperator(m, idim, odim)
end

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T}}, Φ::KrausOperators{T}) where T<:AbstractMatrix{<:Number}
    # TODO : Check if idim or odim
    Stinespring{T}(sum(k⊗ket(i-1, Φ.idim) for (i, k) in enumerate(Φ.matrices)), Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into dynamical matrix.
"""
function Base.convert(::Type{DynamicalMatrix{T}}, Φ::KrausOperators{T}) where T<:AbstractMatrix{<:Number}
    DynamicalMatrix(sum(res(k) * res(k)' for k in Φ.matrices), Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into list of Kraus operators.
"""
function Base.convert(::Type{KrausOperators{T}}, Φ::SuperOperator{T}) where T<:AbstractMatrix{<:Number}
    F = eigfact(Hermitian(reshuffle(Φ.matrix)))
    v = [sqrt(val)*unres(F.vectors[:,i], Φ.odim) for (i, val) in enumerate(F.values)]
    KrausOperators{T}(v, Φ.idim, Φ.odim)
end

function Base.convert(::Type{KrausOperators{T}}, Φ::SuperOperator{T}) where T<:AbstractSparseMatrix{<:Number}
    warn("converting to full matrix")
    convert(KrausOperators{T}, full(Φ))
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into dynamical matrix.
"""
function Base.convert(::Type{DynamicalMatrix{T}}, Φ::SuperOperator{T}) where T<:AbstractMatrix{<:Number}
    DynamicalMatrix{T}(reshuffle(Φ.matrix, [Φ.idim Φ.idim; Φ.odim Φ.odim]), Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T}}, Φ::SuperOperator{T}) where T<:AbstractMatrix{<:Number}
    convert(Stinespring{T},convert(KrausOperators{T}, Φ))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into list of Kraus operators.
"""
function Base.convert(::Type{KrausOperators{T}}, Φ::DynamicalMatrix{T}) where T<:AbstractMatrix{<:Number}
    F = eigfact(Hermitian(Φ))
    v = T[]
    for i in 1:length(F.values)
        if F.values[i] >= 0.0
            push!(v, sqrt(F.values[i]) * unres(F.vectors[:,i], Φ.odim))
        else
            push!(v, zero(unres(F.vectors[:,i], Φ.odim)))
        end
    end
    KrausOperators(v, Φ.idim, Φ.odim)
end

function Base.convert(::Type{KrausOperators{T}}, Φ::DynamicalMatrix{T}) where T<:AbstractSparseMatrix{<:Number}
    warn("converting to full matrix")
    convert(KrausOperators{T}, full(Φ))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T}}, Φ::DynamicalMatrix{T}) where T<:AbstractMatrix{<:Number}
    convert(Stinespring{T}, convert(KrausOperators{T}, Φ))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into super-operator matrix.
"""
function Base.convert(::Type{SuperOperator{T}}, Φ::DynamicalMatrix{T}) where T<:AbstractMatrix{<:Number}
    SuperOperator(reshuffle(Φ.matrix, [Φ.idim Φ.idim; Φ.odim Φ.odim]), Φ.idim, Φ.odim)
end

#### Application of channels
"""
$(SIGNATURES)
- `R`: dynamical matrix.
- `ρ`: quantum state.

Application of dynamical matrix into state `ρ`.
"""
function applychannel(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    unres(reshuffle(Φ.matrix) * res(ρ))
end

"""
$(SIGNATURES)
- `R`: list of Kraus operators.
- `ρ`: quantum state.

Application of list of Kraus operators into state `ρ`.
"""
function applychannel(Φ::KrausOperators{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    sum(k * ρ * k' for k in Φ.matrices)
end

"""
$(SIGNATURES)
- `M`: super-operator matrix.
- `ρ`: quantum state.

Application of super-operator matrix into state `ρ`.
"""
function applychannel(Φ::SuperOperator{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    unres(Φ.matrix * res(ρ))
end

"""
$(SIGNATURES)
- `A`: Stinespring representation of quantum channel.
- `ρ`: quantum state.
- `dims`: dimensions of registers of `ρ`.

Application of Stinespring representation of quantum channel into state `ρ`.
"""
function applychannel(Φ::Stinespring{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    # TODO: Check this function carefully
    dims = s.idim^2, s.odim^2
    ptrace(Φ.matrix * ρ * Φ.matrix', dims, 2)
end

applychannel(Φ::IdentityChannel{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number}) = ρ
# TODO: promote type

function applychannel(Φ::AbstractQuantumOperation, ψ::AbstractVector{<:Number})
    applychannel(Φ, proj(ψ))
end

# TODO : Specialise this fucntion for different quantum ops
function kron(Φ1::T, Φ2::T) where {T<:AbstractQuantumOperation{TM}} where {TM<:AbstractMatrix{<:Number}}
    s1 = SuperOperator{TM}(Φ1)
    s2 = SuperOperator{TM}(Φ2)
    T(SuperOperator{TM}(s1.matrix ⊗ s2.matrix, s1.idim * s2.idim, s1.odim * s2.odim))
end

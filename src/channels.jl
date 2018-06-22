import Base.convert
import Base.size
import Base.length

struct KrausOperators{T<:AbstractMatrix{<:Number}}
    matrices::Vector{T}
    idim::Int
    odim::Int
    function KrausOperators(kl::Vector{T}) where T<:AbstractMatrix{<:Number}
        sizes = [size(k) for k in kl]
        for s in sizes[2:end]
            if s!=sizes[1]
                throw(ArgumentError("Kraus operators list contains matrices of different dimmension"))
            end
        end
        # TODO: check if idims and odims agree with matrix size
        idim, odim = sizes[1]
        new{T}(kl, idim, odim)
    end
end

length(kops::KrausOperators) = length(kops.matrices)
size(kops::KrausOperators) = (idim, odim)

# TODO: create iterator over KrausOperators, length and size function
struct SuperOperator{T<:AbstractMatrix{<:Number}}
    matrix::T
    idim::Int
    odim::Int
    # TODO: write inner constructor
end

struct DynamicalMatrix{T<:AbstractMatrix{<:Number}}
    matrix::T
    idim::Int
    odim::Int
end

struct Stinespring{T<:AbstractMatrix{<:Number}}
    matrix::T
    idim::Int
    odim::Int
end

for _type in (:SuperOperator, :DynamicalMatrix, :Stinespring)
    @eval size(m::$_type) = size(m.matrix)
end

for _type in (:SuperOperator, :DynamicalMatrix, :Stinespring)
    @eval convert(::Type{<:AbstractMatrix{<:Number}}, m::$_type) = m.matrix
end

#### Relationship among representations of channels

"""_
$(SIGNATURES)
- `kraus`: list of Kraus operators.
- `atol`: tolerance of approximation.

Checks if set of Kraus operators fulfill completness relation.
"""
function iscptp(krausops::KrausOperators{<:AbstractMatrix{<:Number}}, atol=1e-08)
    complentess_relation = sum(k'*k for k in krausops.matrices)
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

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into super-operator matrix.
"""
function Base.convert(::Type{SuperOperator{<:AbstractMatrix{<:Number}}}, kops::KrausOperators{<:AbstractMatrix{<:Number}})
    SuperOperator(sum(k⊗(conj.(k)) for k in kops.matrices), kops.idim, kops.odim)
end

"""
$(SIGNATURES)
- `channel`: quantum channel map.
- `dim`: square root of the super-operator matrix dimension.

Transforms quntum channel into super-operator matrix.
"""
function SuperOperator{T}(channel::Function, idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
    odim > 0 ? () : error("Channel dimension has to be nonnegative")

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
function Base.convert(::Type{Stinespring{T}}, kops::KrausOperators{T}) where T<:AbstractMatrix{<:Number}
    # TODO : Check if idim or odim
    Stinespring{T}(sum(k⊗ket(i-1, kops.idim) for (i, k) in enumerate(kops.matrices)), kops.idim, kops.odim)
end

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into dynamical matrix.
"""
function Base.convert(::Type{DynamicalMatrix{T}}, kops::KrausOperators{T}) where T<:AbstractMatrix{<:Number}
    DynamicalMatrix(sum(res(k) * res(k)' for k in kops.matrices), kops.idim, kops.odim)
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into list of Kraus operators.
"""
function Base.convert(::Type{KrausOperators{T}}, sop::SuperOperator{T}) where T<:AbstractMatrix{<:Number}
    F = eigfact(Hermitian(reshuffle(sop.matrix)))
    kops = [sqrt(val)*unres(F.vectors[:,i], sop.odim) for (i, val) in enumerate(F.values)]
    KrausOperators{T}(kops, sop.idim, sop.odim)
end

function Base.convert(::Type{KrausOperators{T}}, sop::SuperOperator{T}) where T<:AbstractSparseMatrix{<:Number}
    warn("converting to full matrix")
    convert(KrausOperators{T}, full(sop))
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into dynamical matrix.
"""
function Base.convert(::Type{DynamicalMatrix{T}}, sop::SuperOperator{T}) where T<:AbstractMatrix{<:Number}
    DynamicalMatrix{T}(reshuffle(sop, [idim idim; odim odim]), sop.idim, sop.odim)
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T}}, sop::SuperOperator{T}) where T<:AbstractMatrix{<:Number}
    convert(Stinespring{T},convert(KrausOperators{T}, sop))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into list of Kraus operators.
"""
function Base.convert(::Type{KrausOperators{T}}, dmat::DynamicalMatrix{T}) where T<:AbstractMatrix{<:Number}
    F = eigfact(Hermitian(dmat))
    kops = T[]
    for i in 1:length(F.values)
        if F.values[i] >= 0.0
            push!(kops, sqrt(F.values[i]) * unres(F.vectors[:,i], dmat.odim))
        else
            push!(kops, zero(unres(F.vectors[:,i], dmat.odim)))
        end
    end
    KrausOperators(kops, dmat.idim, dmat.odim)
end

function Base.convert(::Type{KrausOperators{T}}, dmat::DynamicalMatrix{T}) where T<:AbstractSparseMatrix{<:Number}
    warn("converting to full matrix")
    convert(KrausOperators{T}, full(dmat))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T}}, dmat::DynamicalMatrix{T}) where T<:AbstractMatrix{<:Number}
    convert(Stinespring{T}, convert(KrausOperators{T}, dmat))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into super-operator matrix.
"""
function Base.convert(::Type{SuperOperator{<:AbstractMatrix{<:Number}}}, dmat::DynamicalMatrix{<:AbstractMatrix{<:Number}})
    SuperOperator(reshuffle(dmat, [dmat.idim dmat.idim; dmat.odim dmat.odim]), dmat.idim, dmat.odim)
end

#### Application of channels
"""
$(SIGNATURES)
- `R`: dynamical matrix.
- `rho`: quantum state.

Application of dynamical matrix into state `rho`.
"""
function applychannel(dmat::DynamicalMatrix{<:AbstractMatrix{<:Number}}, rho::AbstractMatrix{<:Number})
    unres(reshuffle(dmat.matrix) * res(rho))
end

"""
$(SIGNATURES)
- `R`: list of Kraus operators.
- `rho`: quantum state.

Application of list of Kraus operators into state `rho`.
"""
function applychannel(kops::KrausOperators{<:AbstractMatrix{<:Number}}, rho::AbstractMatrix{<:Number})
    sum(k * rho * k' for k in kops.matrices)
end

"""
$(SIGNATURES)
- `M`: super-operator matrix.
- `rho`: quantum state.

Application of super-operator matrix into state `rho`.
"""
function applychannel(sop::SuperOperator{<:AbstractMatrix{<:Number}}, rho::AbstractMatrix{<:Number})
    unres(sop.matrix * res(rho))
end

"""
$(SIGNATURES)
- `A`: Stinespring representation of quantum channel.
- `rho`: quantum state.
- `dims`: dimensions of registers of `rho`.

Application of Stinespring representation of quantum channel into state `rho`.
"""
function applychannel(s::Stinespring{<:AbstractMatrix{<:Number}}, rho::AbstractMatrix{<:Number}, dims)
    ptrace(s.matrix * rho * s.matrix', dims, 2)
end

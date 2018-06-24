################################################################################
# Channels definitions and constructors
################################################################################

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

function KrausOperators{T}(v::Vector{T}, idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
    # TODO: Check if idim and odim are compatible with matrix
    KrausOperators(v)
end

length(Φ::KrausOperators) = length(Φ.matrices)
function orthogonalize(Φ::KrausOperators{T}) where {T<:AbstractMatrix{<:Number}}
    KrausOperators{T}(DynamicalMatrix{T}(Φ))
end

# TODO: create iterator over KrausOperators see: https://docs.julialang.org/en/v0.6.3/manual/interfaces/

struct SuperOperator{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrix::T
    idim::Int
    odim::Int
    function SuperOperator(m::T) where T<:AbstractMatrix{<:Number}
        r, c = size(m)
        sr = isqrt(r)
        sc = isqrt(c)
        if r!=sr^2 || c!=sc^2
            throw(ArgumentError("Superoperator matrix has bad dimensions"))
        end
        idim, odim = sr, sc
        new{T}(m, idim, odim)
    end
end

function SuperOperator{T}(m::T, idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
    # TODO: Check if idim and odim are compatible with matrix
    SuperOperator(m)
end

"""
$(SIGNATURES)
- `channel`: quantum channel map.
- `dim`: square root of the super-operator matrix dimension.

Transforms quntum channel into super-operator matrix.
"""
function SuperOperator{T}(channel::Function, idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
    error("Broken")
    odim > 0 ? () : error("Channel dimension has to be nonnegative") # TODO: fix

    m = zeros(T, idim^2, odim^2)
    for (i, e) in enumerate(base_matrices(idim)) # TODO : base_matrices should be not only square
        m[:, i] = res(channel(e))
    end
    SuperOperator(m, idim, odim)
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

struct UnitaryChannel{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrix::T
    idim::Int
    odim::Int
    function UnitaryChannel(m::T) where T<:AbstractMatrix{<:Number}
        r, c = size(m)
        if r!=c
            throw(ArgumentError("UnitaryChannel matrix has to be square"))
        end
        idim, odim = r, c
        new{T}(m, idim, odim)
    end
end

function UnitaryChannel{T}(m::T, idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
    # TODO: Check if idim and odim are compatible with matrix
    UnitaryChannel(m)
end


struct IdentityChannel{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    idim::Int
    odim::Int
end

################################################################################
# size() function
################################################################################
size(Φ::KrausOperators) = (idim, odim)

for qop in (:SuperOperator, :DynamicalMatrix, :Stinespring, :UnitaryChannel)
    @eval size(Φ::$qop) = size(Φ.matrix)
end

size(Φ::IdentityChannel) = (idim, odim)

################################################################################
# making channels callable
################################################################################
for qop in (:KrausOperators, :SuperOperator, :DynamicalMatrix, :Stinespring, :UnitaryChannel)
    @eval begin
        function (Φ::$qop)(ρ)
            applychannel(Φ, ρ)
        end
    end
end

################################################################################
# conversions functions
################################################################################
for qop in (:SuperOperator, :DynamicalMatrix, :Stinespring, :UnitaryChannel)
    @eval convert(::Type{<:AbstractMatrix{<:Number}}, Φ::$qop) = Φ.matrix
end

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into super-operator matrix.
"""
function Base.convert(::Type{SuperOperator{T1}}, Φ::KrausOperators{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    SuperOperator{T1}(sum(k⊗(conj.(k)) for k in Φ.matrices), Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T1}}, Φ::KrausOperators{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    ko = orthogonalize(Φ)
    Stinespring{T1}(sum(k ⊗ ket(i-1, 2*ko.odim) for (i, k) in enumerate(ko.matrices)), ko.idim, ko.odim)
end

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into dynamical matrix.
"""
function Base.convert(::Type{DynamicalMatrix{T1}}, Φ::KrausOperators{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    DynamicalMatrix{T1}(sum(res(k) * res(k)' for k in Φ.matrices), Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into list of Kraus operators.
"""
function Base.convert(::Type{KrausOperators{T1}}, Φ::SuperOperator{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    convert(KrausOperators{T1}, convert(DynamicalMatrix{T2}, Φ))
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into dynamical matrix.
"""
function Base.convert(::Type{DynamicalMatrix{T1}}, Φ::SuperOperator{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    DynamicalMatrix{T1}(reshuffle(Φ.matrix, [Φ.idim Φ.idim; Φ.odim Φ.odim]), Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T1}}, Φ::SuperOperator{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    convert(Stinespring{T1},convert(KrausOperators{T1}, Φ))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into list of Kraus operators.
"""
function Base.convert(::Type{KrausOperators{T1}}, Φ::DynamicalMatrix{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    F = eigfact(Hermitian(Φ.matrix))
    v = T1[]
    for i in 1:length(F.values)
        if F.values[i] >= 0.0
            push!(v, sqrt(F.values[i]) * unres(F.vectors[:,i], Φ.odim))
        else
            push!(v, zero(unres(F.vectors[:,i], Φ.odim)))
        end
    end
    KrausOperators{T1}(v, Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T1}}, Φ::DynamicalMatrix{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    convert(Stinespring{T1}, convert(KrausOperators{T1}, Φ))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into super-operator matrix.
"""
function Base.convert(::Type{SuperOperator{T1}}, Φ::DynamicalMatrix{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    SuperOperator{T1}(reshuffle(Φ.matrix, [Φ.idim Φ.odim; Φ.idim Φ.odim]), Φ.idim, Φ.odim)
end

function Base.convert(::Type{KrausOperators{T1}}, Φ::UnitaryChannel{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    KrausOperators{T1}(T1[Φ], Φ.idim, Φ.odim)
end

for ch in (:SuperOperator, :DynamicalMatrix, :Stinespring)
    @eval begin
        function Base.convert(::Type{$ch{T1}}, Φ::UnitaryChannel{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
            convert($ch{T}, convert(KrausOperators{T}, Φ))
        end
    end
end

################################################################################
# Application of channels
################################################################################
"""
$(SIGNATURES)
- `R`: dynamical matrix.
- `ρ`: quantum state.

Application of dynamical matrix into state `ρ`.
"""
function applychannel(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    ptrace(Φ.matrix * (eye(Φ.idim)⊗transpose(ρ)), [Φ.idim, Φ.odim], [2])
end

"""
$(SIGNATURES)
- `kraus_list`: list of vectors.
- `ρ`: input matrix.

Return mapping of `kraus_list` on `ρ`. Krauss representation of quantum channel
\$\\Phi\$ is a set \$\\{K_i\\}_{i\\in I}\$ of bounded operators on \$\\mathcal{H}\$
such that \$\\sum_{i\\in I} K_i^\\dagger K_i = \\mathcal{1}\$.
Then \$\\Phi(\\rho)=\\sum_{i\\in I} K_i \\rho K_i^\\dagger\$.
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
    ptrace(Φ.matrix * ρ * Φ.matrix', [Φ.idim, Φ.odim^2], [2])
end

function applychannel(Φ::IdentityChannel{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    # TODO: promote type
    ρ
end

function applychannel(Φ::UnitaryChannel{<:AbstractMatrix{<:Number}}, ρ::AbstractMatrix{<:Number})
    # TODO: promote type
    Φ*ρ*Φ'
end

function applychannel(Φ::AbstractQuantumOperation, ψ::AbstractVector{<:Number})
    # TODO: promote type
    applychannel(Φ, proj(ψ))
end

function applychannel(Φ::UnitaryChannel{<:AbstractMatrix{<:Number}}, ψ::AbstractVector{<:Number})
    # TODO: promote type
    Φ*ψ
end

function applychannel(Φ::IdentityChannel{<:AbstractMatrix{<:Number}}, ψ::AbstractVector{<:Number})
    # TODO: promote type
    ψ
end

################################################################################
# composition functions
################################################################################

# TODO : Specialise this function for different quantum ops
function Base.kron(Φ1::T, Φ2::T) where {T<:AbstractQuantumOperation{TM}} where {TM<:AbstractMatrix{<:Number}}
    s1 = SuperOperator{TM}(Φ1)
    s2 = SuperOperator{TM}(Φ2)
    T(SuperOperator{TM}(s1.matrix ⊗ s2.matrix, s1.idim * s2.idim, s1.odim * s2.odim))
end

function Base.kron(Φ1::T, Φ2::T) where {T<:UnitaryChannel{TM}} where {TM<:AbstractMatrix{<:Number}}
    UnitaryChannel(Φ1.matrix ⊗ Φ2.matrix, Φ1.idim * Φ2.idim, Φ1.odim * Φ2.odim)
end

# """
#     Compose channels in sequence
# """
# TODO: allow for composition of different quantum ops - create type hierarchy
function compose(Φ1::T, Φ2::T) where {T<:AbstractQuantumOperation{TM}} where {TM<:AbstractMatrix{<:Number}}
    # TODO : Specialise this function for different quantum ops
    if Φ1.odim != Φ2.idim
        throw(ArgumentError("Channels are incompatible"))
    end
    s1 = SuperOperator{TM}(Φ1)
    s2 = SuperOperator{TM}(Φ2)
    T(SuperOperator{TM}(s1.matrix * s2.matrix, s1.idim, s2.odim))
end

function Base.:*(Φ1::T, Φ2::T) where {T<:AbstractQuantumOperation{TM}} where {TM<:AbstractMatrix{<:Number}}
    compose(Φ1, Φ2)
end

################################################################################
# Channels permutations
################################################################################
# function permutesystems(Φ::T, idims::Vector{Int}, odims::Vector{Int}, perm::Vector{Int}) where {T<:AbstractQuantumOperation{TM}, TM<:AbstractMatrix{<:Number}}
#     error("Not implemented")
#     ko = KrausOperators(Φ)
#     for k in ko
#         #permutesystems(k, sdims, )
#     end
# end

################################################################################
# CPTP, CPTNI
################################################################################
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

function iscptp(Φ::SuperOperator)
    #TODO: implement
end

function iscptp(Φ::DynamicalMatrix)
    #TODO: implement
end

function iscptp(Φ::Stinespring)
    #TODO: implement
end

function iscptp(Φ::UnitaryChannel)
    #TODO: implement
end

function iscp(Φ::SuperOperator)
    #TODO: implement
end

function iscptni(Φ::SuperOperator)
    #TODO: implement
end

function istp(Φ::SuperOperator)
    #TODO: implement
end

function istni(Φ::SuperOperator)
    #TODO: implement
end

export AbstractQuantumOperation, KrausOperators, SuperOperator, DynamicalMatrix,
    Stinespring, UnitaryChannel, IdentityChannel, POVMMeasurement,
    PostSelectionMeasurement, ispovm, iseffect, iscp, istp, istni, iscptp,
    iscptni, applychannel, compose, isidentity, ispositive, represent

################################################################################
# Channels definitions and constructors
################################################################################

abstract type AbstractQuantumOperation{T<:AbstractMatrix{<:Number}} end

"""
$(SIGNATURES)
- `T`: quantum channel map.

Representation of quantum channel by Kraus operators.
"""
struct KrausOperators{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrices::Vector{T}
    idim::Int
    odim::Int
    function KrausOperators{T1}(v::Vector{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
        sizes = [size(k) for k in v]
        for s in sizes[2:end]
            if s!=sizes[1]
                throw(ArgumentError("Kraus operators list contains matrices of different dimension"))
            end
        end
        odim, idim = sizes[1]
        new{T1}(map(T1, v), idim, odim)
    end
end

function KrausOperators{T1}(v::Vector{T2}, idim::Int, odim::Int) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    all((odim, idim) == size(k) for k in v) ? () : throw(ArgumentError("Matrix size and operator dimensions mismatch"))
    KrausOperators{T1}(v)
end

length(Φ::KrausOperators) = length(Φ.matrices)
function orthogonalize(Φ::KrausOperators{T}) where {T<:AbstractMatrix{<:Number}}
    convert(KrausOperators{T}, convert(DynamicalMatrix{T}, Φ))
end

# TODO: create iterator over KrausOperators see: https://docs.julialang.org/en/v0.6.3/manual/interfaces/

"""
$(SIGNATURES)
- `T`: quantum channel map.

Representation of quantum channel by super-operator.
"""
struct SuperOperator{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrix::T
    idim::Int
    odim::Int
    function SuperOperator{T1}(m::T2) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
        r, c = size(m)
        sr = isqrt(r)
        sc = isqrt(c)
        if r!=sr^2 || c!=sc^2
            throw(ArgumentError("Superoperator matrix has invalid dimensions"))
        end
        odim, idim = sr, sc
        new{T1}(convert(T1, m), idim, odim)
    end
end

function SuperOperator{T1}(m::T2, idim::Int, odim::Int) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    (odim^2, idim^2) == size(m) ? () : throw(ArgumentError("Matrix size and operator dimensions mismatch"))
    SuperOperator{T1}(m)
end

"""
$(SIGNATURES)
- `channel`: quantum channel map.
- `idim`: square root of the [super-operator](https://en.wikipedia.org/wiki/Superoperator) matrix input dimension.
- `odim`: square root of the [super-operator](https://en.wikipedia.org/wiki/Superoperator) matrix output dimension.

Transforms quntum channel into super-operator matrix.
"""
function SuperOperator{T}(channel::Function, idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
    odim > 0 && idim > 0 ? () : throw(ArgumentError("Channel dimensions have to be nonnegative"))

    m = zeros(eltype(T), idim^2, odim^2)
    for (i, e) in enumerate(ElementaryBasisIterator{Matrix{Int}}(idim, odim))
        m[:, i] = res(channel(e))
    end
    SuperOperator(m, idim, odim)
end

"""
$(SIGNATURES)
- `T`: quantum channel map.

Representation of quantum channel by Dynamical matrix operators.
"""
struct DynamicalMatrix{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrix::T
    idim::Int
    odim::Int
    function DynamicalMatrix{T1}(m, idim, odim) where {T1<:AbstractMatrix{<:Number}}
        r, c = size(m)
        if r!=c || r!=idim*odim
            throw(ArgumentError("DynamicalMatrix matrix has invalid dimensions"))
        end
        new(convert(T1, m), idim, odim)
    end
end

"""
$(SIGNATURES)
- `T`: quantum channel map.

Stinespring representation of quantum channel.
"""
struct Stinespring{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrix::T
    idim::Int
    odim::Int
    function Stinespring{T1}(m, idim, odim) where {T1<:AbstractMatrix{<:Number}}
        r, c = size(m)
        if r!=idim * (odim^2) || c!=idim
            throw(ArgumentError("Stinespring matrix has invalid dimensions"))
        end
        new(T1(m), idim, odim)
    end
end

"""
$(SIGNATURES)
- `T`: quantum channel map.

Representation of unitary channel.
"""
struct UnitaryChannel{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrix::T
    idim::Int
    odim::Int
    function UnitaryChannel{T1}(m::T2) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
        odim, idim = size(m)
        idim == odim ? () : throw(ArgumentError("UnitaryChannel matrix has to be square"))
        new{T1}(convert(T1,m), idim, odim)
    end
end

function UnitaryChannel{T1}(m::T2, idim::Int, odim::Int) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    (odim, idim) == size(m) ? () : throw(ArgumentError("Matrix size and operator dimensions mismatch"))
    UnitaryChannel{T1}(convert(T1,m))
end

"""
$(SIGNATURES)
- `T`: quantum channel map.

Representation of identity channel.
"""
struct IdentityChannel{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    idim::Int
    odim::Int
    function IdentityChannel{T}(dim::Int) where T<:AbstractMatrix{<:Number}
         new{T}(dim, dim)
    end
end

IdentityChannel(dim::Int) = IdentityChannel{Matrix{ComplexF64}}(dim)

################################################################################
# measurements
################################################################################
struct POVMMeasurement{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrices::Vector{T}
    idim::Int
    odim::Int
    function POVMMeasurement{T1}(v::Vector{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
        sizes = [size(p) for p in v]
        for s in sizes[2:end]
            if s!=sizes[1]
                throw(ArgumentError("POVM operators list contains matrices of different dimension"))
            end
        end
        idim = size(v[1], 1) # TODO or size(v[1], 2) ?
        odim = length(v)

        new{T1}(map(T1, v), idim, odim)
    end
end

function POVMMeasurement{T1}(v::Vector{T2}, idim::Int, odim::Int) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    all((idim, idim) == size(p) for p in v) ? () : throw(ArgumentError("POVMs must be square matrices of size equal to operator inupt dimension"))
    odim == length(v) ? () : throw(ArgumentError("Operator output dimension must match number of POVM operators"))
    POVMMeasurement{T1}(v)
end

function ispovm(Φ::POVMMeasurement{<:AbstractMatrix{<:Number}})
    isidentity(sum(Φ.matrices))
end

struct PostSelectionMeasurement{T<:AbstractMatrix{<:Number}} <: AbstractQuantumOperation{T}
    matrix::T
    idim::Int
    odim::Int
    function PostSelectionMeasurement{T1}(m::T2) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
        odim, idim = size(m)
        new{T1}(convert(T1,m), idim, odim)
    end
end

function PostSelectionMeasurement{T1}(m::T2, idim::Int, odim::Int) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    odim, idim == size(m) ? () : throw(ArgumentError("Matrix size and operator dimensions mismatch"))
    PostSelectionMeasurement{T1}(m)
end

function iseffect(Φ::PostSelectionMeasurement{<:AbstractMatrix{<:Number}})
    e = Φ.matrix
    m = e'*e
    ispositive(one(m) - m)
end

################################################################################
# typeless constructors
################################################################################
for qop in (:SuperOperator, :UnitaryChannel, :PostSelectionMeasurement)
    @eval begin
        function $qop(m::T) where T<:AbstractMatrix{<:Number}
            $qop{T}(m)
        end

        function $qop(m::T, idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
            $qop{T}(m, idim, odim)
        end
    end
end

for qop in (:DynamicalMatrix, :Stinespring)
    @eval begin
        function $qop(m::T, idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
            $qop{T}(m, idim, odim)
        end
    end
end

for qop in (:KrausOperators, :POVMMeasurement)
    @eval begin
        function $qop(v::T) where T<:Vector{M} where M<:AbstractMatrix{<:Number}
            $qop{M}(v)
        end

        function $qop(v::T, idim::Int, odim::Int) where  T<:Vector{M} where M<:AbstractMatrix{<:Number}
            $qop{M}(v)
        end
    end
end


################################################################################
# size() function
################################################################################
# for qop in (:KrausOperators, :POVMMeasurement)
#     @eval size(Φ::$qop) = (Φ.idim, Φ.odim)
# end
#
# for qop in (:SuperOperator, :DynamicalMatrix, :Stinespring,
#             :UnitaryChannel, :PostSelectionMeasurement)
#     @eval size(Φ::$qop) = size(Φ.matrix)
# end
#
# size(Φ::IdentityChannel) = (Φ.idim, Φ.odim)

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

################################################################################
# represent() function
################################################################################
for qop in (:KrausOperators, :POVMMeasurement)
    @eval represent(Φ::$qop) = Φ.matrices
end

for qop in (:SuperOperator, :DynamicalMatrix, :Stinespring,
            :UnitaryChannel, :PostSelectionMeasurement)
    @eval represent(Φ::$qop) = Φ.matrix
end

represent(Φ::IdentityChannel{T}) where T<:Matrix{<:Number} = Matrix{T}(I, Φ.idim, Φ.idim)
represent(Φ::IdentityChannel) = represent(IdentityChannel{Matrix{ComplexF64}}())

################################################################################
# conversions functions
################################################################################

for qop in (:SuperOperator, :DynamicalMatrix, :Stinespring, :UnitaryChannel,
            :PostSelectionMeasurement)
    @eval convert(::Type{<:AbstractMatrix{<:Number}}, Φ::$qop) = Φ.matrix
end

"""
$(SIGNATURES)
- ?: type.
- `Φ`: list of Kraus operators.

Transforms list of Kraus operators into super-operator matrix.
"""
function Base.convert(::Type{SuperOperator{T1}}, Φ::KrausOperators{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    m = sum(k⊗(conj.(k)) for k in Φ.matrices)
    SuperOperator{T1}(convert(T1,m), Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- ?: type.
- `Φ`: list of Kraus operators.

Transforms list of Kraus operators into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T1}}, Φ::KrausOperators{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    ko = orthogonalize(Φ)
    # TODO: improvement: transform to stacking
    m = sum(k ⊗ ket(i, ko.idim*ko.odim) for (i, k) in enumerate(ko.matrices))
    Stinespring{T1}(convert(T1,m), ko.idim, ko.odim)
end

"""
$(SIGNATURES)
- ?: type.
- `Φ`: list of Kraus operators.

Transforms list of Kraus operators into dynamical matrix.
"""
function Base.convert(::Type{DynamicalMatrix{T1}}, Φ::KrausOperators{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    m = sum(res(k) * res(k)' for k in Φ.matrices)
    DynamicalMatrix{T1}(convert(T1,m), Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- ?: type.
- `Φ`: super-operator matrix.

Transforms super-operator matrix into list of Kraus operators.
"""
function Base.convert(::Type{KrausOperators{T1}}, Φ::SuperOperator{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    convert(KrausOperators{T1}, convert(DynamicalMatrix{T2}, Φ))
end

"""
$(SIGNATURES)
- ?: type.
- `Φ`: super-operator matrix.

Transforms super-operator matrix into dynamical matrix.
"""
function Base.convert(::Type{DynamicalMatrix{T1}}, Φ::SuperOperator{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    m = reshuffle(Φ.matrix, [Φ.odim Φ.odim; Φ.idim Φ.idim])
    DynamicalMatrix{T1}(convert(T1,m), Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- ?: type.
- `Φ`: super-operator matrix.

Transforms super-operator matrix into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T1}}, Φ::SuperOperator{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    convert(Stinespring{T1}, convert(KrausOperators{T1}, Φ))
end

"""
$(SIGNATURES)
- ?: type.
- `Φ`: dynamical matrix.

Transforms dynamical matrix into list of Kraus operators.
"""
function Base.convert(::Type{KrausOperators{T1}}, Φ::DynamicalMatrix{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    isnumbernotint(eltype(T1)) ? () : throw(ArgumentError("Kraus operators element type must be subtype of Real or Complex"))

    d = copy(Φ.matrix)
    for i=1:size(d, 1)
        d[i, i] = real(d[i, i])
    end
    F = eigen(Hermitian(d))

    v = T1[]
    for i in 1:length(F.values)
        if F.values[i] >= 0.0
            push!(v, sqrt(F.values[i]) * unres(F.vectors[:,i], Φ.idim))
        else
            push!(v, zero(unres(F.vectors[:,i], Φ.idim)))
        end
    end
    KrausOperators{T1}(v, Φ.idim, Φ.odim)
end

"""
$(SIGNATURES)
- ?: type.
- `Φ`: dynamical matrix.

Transforms dynamical matrix into Stinespring representation of quantum channel.
"""
function Base.convert(::Type{Stinespring{T1}}, Φ::DynamicalMatrix{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    convert(Stinespring{T1}, convert(KrausOperators{T1}, Φ))
end

"""
$(SIGNATURES)
- ?: type.
- `Φ`: dynamical matrix.

Transforms dynamical matrix into super-operator matrix.
"""
function Base.convert(::Type{SuperOperator{T1}}, Φ::DynamicalMatrix{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    m = reshuffle(Φ.matrix, [Φ.odim Φ.idim; Φ.odim Φ.idim])
    SuperOperator{T1}(convert(T1,m), Φ.idim, Φ.odim)
end

function Base.convert(::Type{KrausOperators{T1}}, Φ::UnitaryChannel{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    KrausOperators{T1}(T1[convert(T1,Φ.matrix)], Φ.idim, Φ.odim)
end

function Base.convert(::Type{KrausOperators{T1}}, Φ::IdentityChannel{T2}) where {T1<:AbstractMatrix{N1}, T2<:AbstractMatrix{N2}} where {N1<:Number, N2<:Number}
    N = promote_type(N1, N2)
    # KrausOperators{T1}(T1[eye(N, Φ.idim)], Φ.idim, Φ.odim)
    KrausOperators{T1}(T1[Matrix{N}(I, Φ.idim, Φ.idim)], Φ.idim, Φ.odim)
end

function Base.convert(::Type{KrausOperators{T1}}, Φ::POVMMeasurement{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    # TODO : Verify!
    # v = T1[ket(i-1, Φ.odim)*bra(j-1, Φ.idim)*sqrt(p) for (i, p) in enumerate(Φ.matrices) for j in 1:Φ.idim]
    isnumbernotint(eltype(T1)) ? () : throw(ArgumentError("Kraus operators element type must be subtype of Real or Complex"))
    v = T1[]
    for (i, p) in enumerate(Φ.matrices)
        sqrtp = sqrt(p)
        k = ket(i, Φ.odim)*sum(bra(j, Φ.idim)*sqrtp for j in 1:Φ.idim)
        push!(v, convert(T1,k))
    end
    KrausOperators{T1}(v, Φ.idim, Φ.odim)
end

function Base.convert(::Type{KrausOperators{T1}}, Φ::PostSelectionMeasurement{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    m = Φ.matrix
    v = T1[convert(T1,m)]
    KrausOperators{T1}(v, Φ.idim, Φ.odim)
end

for chout in (:SuperOperator, :DynamicalMatrix, :Stinespring)
    for chin in (:UnitaryChannel, :IdentityChannel, :POVMMeasurement, :PostSelectionMeasurement)
        @eval begin
            function Base.convert(::Type{$chout{T1}}, Φ::$chin{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
                convert($chout{T1}, convert(KrausOperators{T1}, Φ))
            end
        end
    end
end


# function Base.convert(::Type{KrausOperators{T1}}, Φ::SuperOperator{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
#     convert(KrausOperators{T1}, convert(DynamicalMatrix{T2}, Φ))
# end

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
    Φ*ψ
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
# composition functions
################################################################################

# TODO : Specialise this function for different quantum ops
function Base.kron(::Type{T}, Φ1::T1, Φ2::T2) where {T<:AbstractQuantumOperation{M}, T1<:AbstractQuantumOperation{M1}, T2<:AbstractQuantumOperation{M2}} where {M<:AbstractMatrix{<:Number}, M1<:AbstractMatrix{<:Number}, M2<:AbstractMatrix{<:Number}}
    ko1 = convert(KrausOperators{M1}, Φ1)
    ko2 = convert(KrausOperators{M2}, Φ2)
    # M = promote_type(M1, M2)
    v = M[k1⊗k2 for k1 in ko1.matrices for k2 in ko2.matrices]
    ko = KrausOperators{M}(v, Φ1.idim * Φ2.idim, Φ1.odim * Φ2.odim)
    convert(T, ko)
end

function Base.kron(Φ1::T1, Φ2::T2) where {T1<:AbstractQuantumOperation{M1}, T2<:AbstractQuantumOperation{M2}} where {M1<:AbstractMatrix{<:Number}, M2<:AbstractMatrix{<:Number}}
    M = promote_type(M1, M2)
    kron(KrausOperators{M}, Φ1, Φ2)
end

function Base.kron(::Type{T}, Φ1::UnitaryChannel, Φ2::UnitaryChannel) where {T<:AbstractQuantumOperation{M}} where  {M<:AbstractMatrix{<:Number}}
    matrix = convert(M, Φ1.matrix ⊗ Φ2.matrix)
    uc = UnitaryChannel(matrix, Φ1.idim * Φ2.idim, Φ1.odim * Φ2.odim)
    convert(T, uc)
end

function Base.kron(Φ1::UnitaryChannel{M1}, Φ2::UnitaryChannel{M2}) where {M1<:AbstractMatrix{<:Number}, M2<:AbstractMatrix{<:Number}}
    M = promote_type(M1, M2)
    kron(UnitaryChannel{M}, Φ1, Φ2)
end

function Base.kron(::Type{T}, Φ1::IdentityChannel, Φ2::UnitaryChannel) where {T<:AbstractQuantumOperation{M}} where  {M<:AbstractMatrix{<:Number}}
    matrix = convert(M, Matrix(I, Φ1.odim, Φ1.idim) ⊗ Φ2.matrix)
    uc = UnitaryChannel(matrix, Φ1.idim * Φ2.idim, Φ1.odim * Φ2.odim)
    convert(T, uc)
end

function Base.kron(Φ1::IdentityChannel{M1}, Φ2::UnitaryChannel{M2}) where {M1<:AbstractMatrix{<:Number}, M2<:AbstractMatrix{<:Number}}
    M = promote_type(M1, M2)
    kron(UnitaryChannel{M}, Φ1, Φ2)
end

function Base.kron(::Type{T}, Φ1::UnitaryChannel, Φ2::IdentityChannel) where {T<:AbstractQuantumOperation{M}} where  {M<:AbstractMatrix{<:Number}}
    matrix = convert(M, Φ1.matrix ⊗ Matrix(I, Φ2.odim, Φ2.idim))
    uc = UnitaryChannel(matrix, Φ1.idim * Φ2.idim, Φ1.odim * Φ2.odim)
    convert(T, uc)
end

function Base.kron(Φ1::UnitaryChannel{M1}, Φ2::IdentityChannel{M2}) where {M1<:AbstractMatrix{<:Number}, M2<:AbstractMatrix{<:Number}}
    M = promote_type(M1, M2)
    kron(UnitaryChannel{M}, Φ1, Φ2)
end


#
# function ⊗(Φ1::T1, Φ2::T2) where {T1<:AbstractQuantumOperation{M1}, T2<:AbstractQuantumOperation{M2}} where {M1<:AbstractMatrix{<:Number}, M2<:AbstractMatrix{<:Number}}
#     M = promote_type(M1, M2)
#     kron(SuperOperator{M}, ϕ1, ϕ2)
# end
#
# const ⊗ᵪ = kron_chi


# # """
# #     Compose channels in sequence
# # """
# # TODO: allow for composition of different quantum ops - create type hierarchy
# function compose(Φ1::T, Φ2::T) where {T<:AbstractQuantumOperation{TM}} where {TM<:AbstractMatrix{<:Number}}
#     # TODO : Specialise this function for different quantum ops
#     if Φ1.odim != Φ2.idim
#         throw(ArgumentError("Channels are incompatible"))
#     end
#     s1 = SuperOperator{TM}(Φ1)
#     s2 = SuperOperator{TM}(Φ2)
#     T(SuperOperator{TM}(s1.matrix * s2.matrix, s1.idim, s2.odim))
# end

function compose(::Type{T}, Φ1::T1, Φ2::T2) where {T<:AbstractQuantumOperation{M}, T1<:AbstractQuantumOperation{M1}, T2<:AbstractQuantumOperation{M2}} where {M<:AbstractMatrix{<:Number}, M1<:AbstractMatrix{<:Number}, M2<:AbstractMatrix{<:Number}}
    if Φ1.odim != Φ2.idim
        throw(ArgumentError("Channels are incompatible"))
    end
    s1 = convert(SuperOperator{M1}, Φ1)
    s2 = convert(SuperOperator{M2}, Φ2)
    matrix = convert(M, s1.matrix * s2.matrix)
    so = SuperOperator{M}(matrix, s1.idim, s2.odim)
    convert(T, so)
end

function compose(Φ1::T, Φ2::T) where {T<:AbstractQuantumOperation{<:Number}}
    compose(T, Φ1, Φ2)
end

function Base.:*(Φ1::T, Φ2::T) where {T<:AbstractQuantumOperation{<:Number}}
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
# Channels printing
################################################################################
function Base.show(io::IO, Φ::AbstractQuantumOperation{<:Matrix{<:Number}})
    println(io, typeof(Φ))
    println(io, "    dimensions: ($(Φ.idim), $(Φ.odim))")
    if :matrix in fieldnames(typeof(Φ))
        print(io, "    ")
        print(io, Φ.matrix)
    end
    if :matrices in fieldnames(typeof(Φ))
        for (i,m) in enumerate(Φ.matrices)
            print(io, "    ")
            print(io, m)
            if i < length(Φ.matrices) println(io, "") end
        end
    end
end
# Base.show(io::IO, m::MIME"text/plain", Φ::AbstractQuantumOperation{<:Matrix{<:Number}}) = show(io, m, Φ)

################################################################################
# CP
################################################################################
"""
$(SIGNATURES)
- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is completely positive.
"""
function iscp end

function iscp(Φ::KrausOperators{<:AbstractMatrix{<:Number}}; atol=1e-13)
    # by definition Kraus operators represent a CP map
    true
end

function iscp(Φ::SuperOperator{T}; atol=1e-13) where T<:AbstractMatrix{<:Number}
    iscp(convert(DynamicalMatrix{T}, Φ), atol=atol)
end

function iscp(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}; atol=1e-13)
    ispositive(Φ.matrix)
end

function iscp(Φ::Stinespring{<:AbstractMatrix{<:Number}}; atol=1e-13)
    # by definition Stinespring operator represents a CP map(?)
    true
end

function iscp(Φ::UnitaryChannel; atol=1e-13)
    # by definition Unitary operator represents a CP map
    true
end


################################################################################
# TNI
################################################################################
"""
$(SIGNATURES)
- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is trace non-increasing.
"""
function istni end

function istni(Φ::KrausOperators{<:AbstractMatrix{<:Number}}; atol=1e-13)
    cr = sum(k'*k for k in Φ.matrices)
    ispositive(one(cr) - cr, atol=atol)
end

function istni(Φ::SuperOperator{T}; atol=1e-13) where T<:AbstractMatrix{<:Number}
    istni(convert(DynamicalMatrix{T}, Φ))
end

function istni(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}; atol=1e-13)
    pt = ptrace(Φ.matrix, [Φ.odim, Φ.idim], [1])
    ispositive(one(pt) - pt, atol=atol)
end

function istni(Φ::Stinespring{<:AbstractMatrix{<:Number}}; atol=1e-13)
    u = Φ.matrix
    m = u'*u
    ispositive(one(m) - m, atol=atol)
end

function istni(Φ::UnitaryChannel; atol=1e-13)
    iscptp(Φ, atol=atol)
end

################################################################################
# TP
################################################################################
"""
$(SIGNATURES)
- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is trace preserving.
"""
function istp end

function istp(Φ::KrausOperators{<:AbstractMatrix{<:Number}}; atol=1e-13)
    cr = sum(k'*k for k in Φ.matrices)
    isidentity(cr, atol=atol)
end

function istp(Φ::SuperOperator{T}; atol=1e-13) where T<:AbstractMatrix{<:Number}
    istp(convert(DynamicalMatrix{T}, Φ), atol=atol)
end

function istp(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}; atol=1e-13)
    pt = ptrace(Φ.matrix, [Φ.odim, Φ.idim], [1])
    isidentity(pt, atol=atol)
end

function istp(Φ::Stinespring{<:AbstractMatrix{<:Number}}; atol=1e-13)
    u = Φ.matrix
    isidentity(u'*u, atol=atol)
end

function istp(Φ::UnitaryChannel; atol=1e-13)
    u = Φ.matrix
    isidentity(u'*u, atol=atol) && isidentity(u*u', atol=atol)
end

################################################################################
# CPTP, CPTNI
################################################################################
function iscptp(Φ::AbstractQuantumOperation; atol=1e-13)
    iscp(Φ, atol=atol) && istp(Φ, atol=atol)
end

function iscptni(Φ::AbstractQuantumOperation; atol=1e-13)
    iscp(Φ, atol=atol) && istni(Φ, atol=atol)
end

export AbstractQuantumOperation, KrausOperators, SuperOperator, DynamicalMatrix,
    Stinespring, UnitaryChannel, IdentityChannel, POVMMeasurement,
    PostSelectionMeasurement
    
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
        idim = size(v[1], 1) 
        odim = length(v)

        new{T1}(map(T1, v), idim, odim)
    end
end

function POVMMeasurement{T1}(v::Vector{T2}, idim::Int, odim::Int) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    all((idim, idim) == size(p) for p in v) ? () : throw(ArgumentError("POVMs must be square matrices of size equal to operator inupt dimension"))
    odim == length(v) ? () : throw(ArgumentError("Operator output dimension must match number of POVM operators"))
    POVMMeasurement{T1}(v)
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
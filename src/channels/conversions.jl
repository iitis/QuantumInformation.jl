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
        push!(v, convert(T1, k))
    end
    KrausOperators{T1}(v, Φ.idim, Φ.odim)
end

function Base.convert(::Type{KrausOperators{T1}}, Φ::PostSelectionMeasurement{T2}) where {T1<:AbstractMatrix{<:Number}, T2<:AbstractMatrix{<:Number}}
    m = Φ.matrix
    v = T1[convert(T1, m)]
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
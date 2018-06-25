function Base.convert(::Type{KrausOperators{T1}}, Φ::SuperOperator{T2}) where {T1<:AbstractSparseMatrix{<:Number}, T2<:AbstractSparseMatrix{<:Number}}
    # FIXME
    # TODO : Try to implement it without conversion
    warn("convertion from sparse to sparse using full matrix")
    fullΦ = full(Φ)
    T = typeof(fullΦ)
    ko = convert(KrausOperators{T}, fullΦ)
    KrausOperators{T1}([sparse(k) for k in ko.matrices], Φ.idim, Φ.odim)
end

# TODO : add sparse to full and full to sparse

function Base.convert(::Type{KrausOperators{T}}, Φ::DynamicalMatrix{T}) where T<:AbstractSparseMatrix{<:Number}
    warn("converting to full matrix")
    convert(KrausOperators{T}, full(Φ))
end

# TODO : add sparse to full and full to sparse

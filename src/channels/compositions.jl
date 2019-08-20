export compose, kron

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
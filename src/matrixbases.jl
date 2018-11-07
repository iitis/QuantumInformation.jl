abstract type AbstractMatrixBase{T<:AbstractMatrix{<:Number}} end
struct ElementaryMatrixBase{T<:AbstractMatrix{<:Number}} <: AbstractMatrixBase{T}
    dim::Int
    length::Int
end

function ElementaryMatrixBase{T}(dim::Int) where T <: AbstractMatrix{<:Number}
    dim > 0 ? () : error("Operator dimension has to be nonnegative")
    ElementaryMatrixBase{T}(dim, dim^2)
end

Base.length(iter::AbstractMatrixBase{T}) where T<:AbstractMatrix{<:Number} = iter.length
Base.eltype(iter::AbstractMatrixBase{T}) where T<:AbstractMatrix{<:Number} = T


# function ElementaryMatrixBase{T}(::Type{T}, dim::Int) where T<:AbstractMatrix{<:Number}
#     ElementaryMatrixBase(base_matrices(T, dim), dim)
# end

# TODO: allow rectangular matrices
function Base.iterate(iter::ElementaryMatrixBase{M}, state=(Iterators.product(1:iter.dim,1:iter.dim),)) where M<:AbstractMatrix{T} where T<:Number
    # for i=1:dim, j=1:dim
    #     push!(c, )
    # end
    it = first(state)
    i, j = it
    element = ketbra(T, j, i, iter.dim)
    if i == iter.dim && j == iter.dim
        return nothing
    end

    state = iterate(it, state)
    return(element, state)
end

import Base.length

import Base.collect


"""
$(SIGNATURES)
- `dim`: length of the matrix.

Returns elementary matrices of dimension `dim` x `dim`.
"""
base_matrices(dim::Int) = base_matrices(ComplexF64, dim)

"""
$(SIGNATURES)
- `dim`: dimensions of registers of `ρ`.

Returns elementary hermitian matrices of dimension dim x dim.
"""
base_hermitian_matrices(dim) = Channel(ctype=Matrix{ComplexF64}) do bhm
    for (a, b) in Base.product(0:dim-1, 0:dim-1)
        if a > b
            x = 1 / sqrt(2) * (1im * ketbra(a, b, dim) - 1im * ketbra(b, a, dim))
            push!(bhm, x)
        elseif a < b
            x = 1 / sqrt(2) * (ketbra(a, b, dim) + ketbra(b, a, dim))
            push!(bhm, x)
        else
            x = ketbra(a, b, dim)
            push!(bhm, x)
        end
    end
    bhm
end

base_generlized_pauli_matrices(d) = Channel(ctype=Matrix{ComplexF64}) do gm
    E(i,j,d) = ketbra(j,i,d)

    sm = sum([E(i,i,d) for i in 0:d-1])
    push!(gm, sm)
    for j in 0:d-2
        for k in j+1:d-1
            sm = E(j,k,d) + E(k,j,d)
            push!(gm, sm)
        end
    end

    for j in 0:d-2
        for k in j+1:d-1
            sm = -1im*(E(j,k,d) - E(k,j,d))
            push!(gm, sm)
        end
    end

    for l in 1:d-1
        sm = sqrt(2.0/(l*(l+1)))*(sum([E(j-1,j-1,d) for j in 1:l]) - l*E(l,l,d))
        push!(gm, sm)
    end
    gm
end

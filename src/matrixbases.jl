import .Base: length, iterate
export AbstractMatrixBasisIterator, HermitianBasisIterator, AbstractBasis,
    AbstractMatrixBasis, HermitianBasis, hermitianbasis, represent, combine
abstract type AbstractMatrixBasisIterator{T<:AbstractMatrix} end
struct HermitianBasisIterator{T} <: AbstractMatrixBasisIterator{T} 
    dim::Int
end

abstract type AbstractBasis end
abstract type AbstractMatrixBasis{T} <: AbstractBasis where T<:AbstractMatrix{<:Number} end 

struct HermitianBasis{T} <: AbstractMatrixBasis{T} 
    iterator::HermitianBasisIterator{T}

    function HermitianBasis{T}(dim::Integer) where T<:AbstractMatrix{<:Number}
        new(HermitianBasisIterator{T}(dim)) 
    end
end

"""
$(SIGNATURES)
- `dim`: dimensions of the matrix.

Returns elementary hermitian matrices of dimension `dim` x `dim`.
"""
hermitianbasis(T::Type{<:AbstractMatrix{<:Number}}, dim::Int) = HermitianBasisIterator{T}(dim)

hermitianbasis(dim::Int) = hermitianbasis(Matrix{ComplexF64}, dim)


function iterate(itr::HermitianBasisIterator{T}, state=(1,1)) where T<:AbstractMatrix{<:Number}
    dim = itr.dim
    (a, b) = state
    a > dim && return nothing

    if a > b
        x = 1 / sqrt(2) * (1im * ketbra(a, b, dim) - 1im * ketbra(b, a, dim))
    elseif a < b
        x = 1 / sqrt(2) * (ketbra(a, b, dim) + ketbra(b, a, dim))
    else
        x = ketbra(a, b, dim)
    end
    return x, b==dim ? (a+1, 1) : (a, b+1)
end

length(itr::HermitianBasisIterator) = itr.dim^2

function represent(basis::T, m::Matrix{<:Number}) where T<:AbstractMatrixBasis
    tr.([m] .* basis.iterator) 
end

function combine(basis::T, v::Vector{<:Number}) where T<:AbstractMatrixBasis
    sum(basis.iterator .* v)
end
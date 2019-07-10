abstract type AbstractMatrixBasisIterator{T<:AbstractMatrix} end
struct HermitianBasisIterator{T} <: AbstractMatrixBasisIterator{T} 
    dim::Int
end
import .Base: length, iterate
export hermitianbasis

hermitianbasis(dim::Int) = hermitianbasis(Matrix{ComplexF64}, dim)

"""
$(SIGNATURES)
- `dim`: dimensions of the matrix.

Returns elementary hermitian matrices of dimension `dim` x `dim`.
"""
hermitianbasis(T::Type{<:AbstractMatrix{<:Number}}, dim::Int) = HermitianBasisIterator{T}(dim)

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

function in_matrix_basis(basis::T, m::Matrix{<:Number}) where T<:AbstractMatrixBasisIterator
    tr.([m] .* basis) 
end

function from_matrix_basis(basis::T, v::Vector{<:Number}) where T<:AbstractMatrixBasisIterator
    sum(basis .* v)
end
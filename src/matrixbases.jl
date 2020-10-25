import .Base: length, iterate

export AbstractMatrixBasisIterator, HermitianBasisIterator, ElementaryBasisIterator, AbstractBasis,
    AbstractMatrixBasis, HermitianBasis, ElementaryBasis,
    ChannelBasisIterator, AbstractChannelBasis, ChannelBasis, 
    hermitianbasis, channelbasis,
    represent, combine

abstract type AbstractMatrixBasisIterator{T<:AbstractMatrix} end
struct HermitianBasisIterator{T} <: AbstractMatrixBasisIterator{T} 
    dim::Int
end

struct ElementaryBasisIterator{T} <: AbstractMatrixBasisIterator{T} 
    idim::Int
    odim::Int
end

abstract type AbstractBasis end
abstract type AbstractMatrixBasis{T} <: AbstractBasis where T<:AbstractMatrix{<:Number} end 

struct HermitianBasis{T} <: AbstractMatrixBasis{T} 
    iterator::HermitianBasisIterator{T}

    function HermitianBasis{T}(dim::Integer) where T<:AbstractMatrix{<:Number}
        new(HermitianBasisIterator{T}(dim)) 
    end
end

struct ElementaryBasis{T} <: AbstractMatrixBasis{T} 
    iterator::ElementaryBasisIterator{T}

    function ElementaryBasis{T}(idim::Integer, odim::Integer) where T<:AbstractMatrix{<:Number}
        new(ElementaryBasisIterator{T}(idim, odim)) 
    end
end

function ElementaryBasis{T}(dim::Integer) where T<:AbstractMatrix{<:Number}
    ElementaryBasisIterator{T}(dim, dim)
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

    Tn = eltype(T)
    if a > b
        x =  (im * ketbra(T, a, b, dim) - im * ketbra(T, b, a, dim)) / sqrt(Tn(2))
    elseif a < b
        x = (ketbra(T, a, b, dim) + ketbra(T, b, a, dim)) / sqrt(Tn(2))
    else
        x = ketbra(T, a, b, dim)
    end
    return x, b==dim ? (a+1, 1) : (a, b+1)
end

length(itr::HermitianBasisIterator) = itr.dim^2

function iterate(itr::ElementaryBasisIterator{T}, state=(1,1)) where T<:AbstractMatrix{<:Number}
    idim, odim = itr.idim, itr.odim
    (a, b) = state
    a > idim  && return nothing

    x = zeros(eltype(T), odim, idim)
    x[b, a] = one(eltype(T))
    return x, b == odim ? (a+1, 1) : (a, b+1)
end

length(itr::ElementaryBasisIterator) = itr.idim * itr.odim

function represent(basis::T, m::Matrix{<:Number}) where T<:AbstractMatrixBasis
    real.(tr.([m] .* basis.iterator))
end

function represent(basis::Type{T}, m::Matrix{<:Number}) where T<:AbstractMatrixBasis
    d = size(m, 1)
    represent(basis{typeof(m)}(d), m)
end

function combine(basis::AbstractMatrixBasis{T}, v::Vector{<:Number}) where T<:AbstractMatrix{<:Number}
    sum(basis.iterator .* v)
end

"""
$(SIGNATURES)

"""
struct ChannelBasisIterator{T} <: AbstractMatrixBasisIterator{T} 
    idim::Int 
    odim::Int
    hitr::HermitianBasisIterator{T}
    function ChannelBasisIterator{T}(idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
        new(idim, odim, HermitianBasisIterator{T}(idim))
    end
end

abstract type AbstractChannelBasis{T} <: AbstractMatrixBasis{T} end 
struct ChannelBasis{T} <: AbstractChannelBasis{T} 
    iterator::ChannelBasisIterator{T}
    function ChannelBasis{T}(idim::Integer, odim::Integer) where T<:AbstractMatrix{<:Number}
        new(ChannelBasisIterator{T}(idim, odim)) 
    end
end
channelbasis(T::Type{<:AbstractMatrix{<:Number}}, idim::Int, odim::Int=idim) = ChannelBasis{T}(idim, odim)

channelbasis(idim::Int, odim::Int=idim) = channelbasis(Matrix{ComplexF64}, idim, odim)

function iterate(itr::ChannelBasisIterator{T}, state=(1,1,1,1)) where T<:AbstractMatrix{<:Number}
    (idim, odim) = (itr.idim, itr.odim)
    hitr = itr.hitr
    (a, c, b, d) = state
    (a == odim && c == odim && d == 2) && return nothing

    Tn = eltype(T)
    if a > c
        x = (ketbra(T, a, c, odim) ⊗ ketbra(T, b, d, idim) + ketbra(T, c, a, odim) ⊗ ketbra(T, d, b, idim)) / sqrt(Tn(2))
    elseif a < c  
        x = (im * ketbra(T, a, c, odim) ⊗ ketbra(T, b, d, idim) - im * ketbra(T, c, a, odim) ⊗ ketbra(T, d, b, idim)) / sqrt(Tn(2))
    elseif a < odim 
        H = iterate(hitr, (b, d))[1]
        x = (diagm(0 => vcat(ones(Tn, a), Tn[-a], zeros(Tn, odim - a-1))) ⊗ H) / sqrt(Tn(a + a^2))
    else 
        x = Matrix{Tn}(I, idim * odim, idim * odim) / sqrt(Tn(idim * odim))
    end
    if d < idim
        newstate = (a, c, b, d+1)
    elseif d == idim && b < idim
        newstate = (a, c, b+1, 1)
    elseif d == idim && b == idim && c < odim
        newstate = (a, c+1, 1, 1)
    else  
        newstate = (a+1, 1, 1, 1)
    end 
    return x, newstate
end 

length(itr::ChannelBasisIterator) = itr.idim^2 * itr.odim^2 - itr.idim^2 + 1

function represent(basis::AbstractChannelBasis{T1}, Φ::AbstractQuantumOperation{T2}) where T1<:AbstractMatrix{<:Number} where T2<:AbstractMatrix{<:Number}
    J = convert(DynamicalMatrix{T1}, Φ)
    represent(basis, J.matrix) 
end

function combine(basis::AbstractChannelBasis{T}, v::Vector{<:Number}) where T<:AbstractMatrix
    m = sum(basis.iterator .* v)
    DynamicalMatrix{T}(m, basis.iterator.idim, basis.iterator.odim)
end
module ElementaryArrays
using LinearAlgebra
import Base: *, +, -, /, size, getindex, setindex!
export ElementaryArray

include("../../src/utils.jl") # FIXME: UGLY HACK



# Array interface https://docs.julialang.org/en/v1/manual/interfaces/index.html#man-interface-array-1
# Linear Algebra https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/index.html
# Permutation groups http://nemocas.github.io/Nemo.jl/v0.6.3/perm/
struct ElementaryArray{T, N} <: AbstractArray{T, N}
    value::T
    idx::NTuple{N, Int}
    dims::NTuple{N, Int}
    function ElementaryArray{T, N}(value::T, idx::NTuple{N, Int}, dims::Vararg{Int}) where {T<:Number, N}
        all(1 <= i <= d for (i, d) in zip(idx, dims)) ? () : throw(BoundsError("index $idx does not fit in dimensions $dims"))
        new(value, idx, dims)
    end
end

function ElementaryArray{T}(idx::NTuple{N, Int}, dims::Vararg{Int}) where {T<:Number, N}
    ElementaryArray(one(T), idx, dims...)
end

function ElementaryArray(value::T, idx::NTuple{N, Int}, dims::Vararg{Int}) where {T<:Number, N}
    ElementaryArray{T, N}(one(T), idx, dims...)
end

const ElementaryVector{T} = ElementaryArray{T, 1}
const ElementaryMatrix{T} = ElementaryArray{T, 2}

# function ElementaryMatrix{T}(value::T, idx::NTuple{2, Int}, dim1, dim2) where T<:Number
#     ElementaryArray{T, 2}(value, idx, dim1, dim2)
# end

# https://docs.julialang.org/en/v1/manual/arrays/index.html#Array-traits-1
# Base.IndexStyle(::Type{<:ElementaryArray}) = IndexLinear()

function Base.getindex(e::ElementaryArray{T, N}, idx::Vararg{Int, N}) where {T<:Number, N} 
    all(1 <= i <= d for (i, d) in zip(idx, e.dims)) ? () : throw(BoundsError("index $idx does not fit in dimensions $(e.dims)"))
    idx == e.idx ? e.value : zero(T)
end

Base.eltype(e::ElementaryArray{T, N}) where {T<:Number, N} = T

Base.size(e::ElementaryArray) = e.dims

function issquare(e::ElementaryMatrix)
    e.dims[1] == e.dims[2] ? true : false
end

function LinearAlgebra.issymmetric(e::ElementaryMatrix)
    issquare(e) ? () : return false
    e.idx[1] == e.idx[2] ? true : false
end

function LinearAlgebra.isposdef(e::ElementaryMatrix)
    issymmetric(e) ? () : return false
    e.value > 0 ? true : false
end

# ! returned type instability
function Base.:+(e1::ElementaryArray{T, N}, e2::ElementaryArray{T, N}) where {T<:Number, N} 
    size(e1) == size(e2) ? () : throw(ArgumentError("ElementaryArrays are of different size"))
    if e1.idx == e2.idx
        value = convert(T, e1.value + e2.value)
        return ElementaryArray{T, N}(value, e1.idx, e1.dims...)
    else
        x = zeros(T, e1.dims)
        x[e1.idx...] += e1.value
        x[e2.idx...] += e2.value
        x  
    end 
end

LinearAlgebra.det(e::ElementaryMatrix) = zero(eltype(e))

LinearAlgebra.logdet(e::ElementaryMatrix) = log(det(e))

function LinearAlgebra.tr(e::ElementaryMatrix)
    issquare(e) ? () : throw(ArgumentError("matrix is not square"))
    e.idx[1] == e.idx[2] ? e.value : zero(eltype(e))
end

function Base.inv(e::ElementaryMatrix)
end

function LinearAlgebra.inv(e::ElementaryMatrix)
end

function Base.adjoint(e::ElementaryMatrix{T}) where T <: Number
    ElementaryMatrix{T}(conj(e.value), reverse(e.idx), reverse(e.dims)...)
end

function Base.transpose(e::ElementaryMatrix{T}) where T <: Number
    ElementaryMatrix{T}(e.value, reverse(e.idx), reverse(e.dims)...)
end

function Base.conj(e::ElementaryArray{N, T}) where {T <: Number, N}
    ElementaryMatrix{T}(conj(e.value), e.idx, e.dims...)
end

Base.Matrix(e::ElementaryMatrix) = Array(e)

function Base.kron(e1::ElementaryMatrix, e2::AbstractMatrix)
end

function Base.kron(e1::AbstractMatrix, e2::ElementaryMatrix)
end

# function index _____

function Base.kron(e1::ElementaryArray{T1, N}, e2::ElementaryArray{T2, N}) where {T1 <: Number, T2 <: Number, N}
    dims = tuple((e1.dims[d] * e2.dims[d] for d in 1:N)...)
    idx = tuple((mixedradix2number((e1.idx[d]-1, e2.idx[d]-1), (e1.dims[d], e2.dims[d])) + 1 for d in 1:N)...)
    
    T = promote_type(T1, T2)
    value = convert(T, e1.value * e2.value)
    ElementaryArray{T, N}(value, idx, dims...)
end

function Base.:*(e1::ElementaryMatrix{T1}, e2::ElementaryMatrix{T2}) where {T1 <: Number, T2 <: Number}
    T = promote_type(T1, T2)
    dims = (e1.dims[1], e2.dims[2])
    value = (e1.idx[2] == e2.idx[1] ? convert(T, e1.value * e2.value) : zero(T))
    idx = (e1.idx[1], e2.idx[2])
    ElementaryMatrix{T}(value, idx, dims...)
end

function Base.reshape(e::ElementaryArray{T, N}, dims::Vararg{Int}) where {T <: Number, N}
    size = prod(e.dims)
    size == prod(dims) ? () : throw(DimensionMismatch("new dimensions $dims must be consistent with array size $size"))
    # TODO look at this monstrosity! it has to be put into a function
    li = mixedradix2number(reverse(tuple((e.idx[d] - 1 for d in 1:N)...)), reverse(e.dims)) 
    idx = tuple((i + 1 for i in number2mixedradix(li, tuple(dims...)))...)
    ElementaryMatrix{T}(e.value, reverse(idx), dims...)
end

LinearAlgebra.lmul!(a::Number, e::ElementaryMatrix) = (e.value *= a)

LinearAlgebra.rmul!(e::ElementaryMatrix, a::Number) = (e.value *= a)

function Base.:*(a::T1, e::ElementaryArray{T2, N}) where {T1<:Number, T2<:Number, N}
    T = promote_type(T1, T2)
    value = convert(T, a * e.value)
    ElementaryArray{T, N}(value, e.idx, e.dims...)
end

Base.:*(e::ElementaryArray, a::Number) = Base.:*(a, e)

# for op in  (:+, :-)
#     @eval begin
#         function $op(n::T1, e::ElementaryArray{T2, N}) where {T1<:Number, T2<:Number, N} 
#             T = promote_type(T1, T2)
#             x = zeros(T1, e.dims)
#             x[e.idx...] $op= n
#             x
#         end
#     end
# end

# for op in  (:+, :-, :*, :/)
#     @eval begin
#         function $op(e::ElementaryArray{T2, N}, n::T1) where {T1<:Number, T2<:Number, N} 
#             T = promote_type(T1, T2)
#             x = zeros(T1, e.dims)
#             x[e.idx...] = $op(1, n)
#             x
#         end
#     end
# end



end


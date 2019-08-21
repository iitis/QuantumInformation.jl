module ElementaryArrays
using LinearAlgebra
import Base: *, +, -, /, size, getindex, setindex!
export ElementaryArray

# Array interface https://docs.julialang.org/en/v1/manual/interfaces/index.html#man-interface-array-1
# Linear Algebra https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/index.html
# Permutation groups http://nemocas.github.io/Nemo.jl/v0.6.3/perm/
struct ElementaryArray{T, N} <: AbstractArray{T, N}
    idx::NTuple{N, Int}
    dims::NTuple{N, Int}
    function ElementaryArray{T, N}(idx::NTuple{N, Int}, dims::Vararg{Int}) where {T<:Number, N}
        all(i <= d for (i, d) in zip(idx, dims)) ? () : throw(ArgumentError("an index is larger that dimmension"))
        new(idx, dims)
    end
end

size(A::ElementaryArray) = A.dims

getindex(A::ElementaryArray{T, N}, I::Vararg{Int, N}) where {T<:Number, N} = I == A.idx ? one(T) : zero(T)

const ElementaryVector{T} = ElementaryArray{T, 1}
const ElementaryMatrix{T} = ElementaryArray{T, 2}

function +(e1::ElementaryArray{T, N}, e2::ElementaryArray{T, N}) where {T<:Number, N} 
    size(e1) == size(e2) ? () : throw(ArgumentError("ElementaryArrays are of different size"))
    x = zeros(T, e1.dims)
    x[e1.idx...] += one(T)
    x[e2.idx...] += one(T)
    x    
end

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

for op in  (:+, :-, :*, :/)
    @eval begin
        function $op(e::ElementaryArray{T2, N}, n::T1) where {T1<:Number, T2<:Number, N} 
            T = promote_type(T1, T2)
            x = zeros(T1, e.dims)
            x[e.idx...] = $op(1, n)
            x
        end
    end
end

@show e1 = ElementaryArray{Float64, 2}((1, 2), 2, 2)
@show e2 = ElementaryArray{Float64, 2}((2, 2), 2, 2)
@show e1 + e2
@show 1 + e2
@show  e2 / 2
end


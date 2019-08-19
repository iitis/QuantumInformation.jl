module ElementaryArrays
using LinearAlgebra
import Base: *, +, -, size, getindex, setindex!
export ElementaryArray

# Array interface https://docs.julialang.org/en/v1/manual/interfaces/index.html#man-interface-array-1
# Linear Algebra https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/index.html
# Permutation groups http://nemocas.github.io/Nemo.jl/v0.6.3/perm/
struct ElementaryArray{T, N} <: AbstractArray{T, N}
    idxs::Set{NTuple{N, Int}}
    dims::NTuple{N, Int}
    function ElementaryArray{T, N}(idxs::Set{NTuple{N, Int}}, dims::Vararg{Int}) where {T<:Number, N}
        @show idxs
        @show dims, typeof(dims)
        # TODO: check dims

        new(idxs, dims)
    end
end

size(A::ElementaryArray) = A.dims

getindex(A::ElementaryArray{T, N}, I::Vararg{Int, N}) where {T<:Number, N} = I in A.idxs ? one(T) : zero(T)

# @show ElementaryArray{Float64, 2}(Set([(1,1),(1,2)]), 2, 2)
end


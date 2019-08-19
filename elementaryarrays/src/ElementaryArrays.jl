module ElementaryArrays
using LinearAlgebra
import Base: *, +, -, size, getindex, setindex!

# Array interface https://docs.julialang.org/en/v1/manual/interfaces/index.html#man-interface-array-1
# Linear Algebra https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/index.html
# Permutation groups http://nemocas.github.io/Nemo.jl/v0.6.3/perm/
struct ElementaryArray{T, N} <: AbstractArray{T, N}
    idxs
    dims
end
end

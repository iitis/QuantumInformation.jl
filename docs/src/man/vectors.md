# Vectors and matrices in Julia

A basic construction of vector in Julia creates a full one-index array containing elements of a number type as presented below.
```julia
julia> x = [0.0, 1.0im]
2-element Array{Complex{Float64},1}:
 0.0+0.0im
 0.0+1.0im
```
A transposition of a column vector return an object of type `RowVector` as shown below
```julia
julia> xt = transpose(x)
1×2 RowVector{Complex{Float64},Array{Complex{Float64},1}}:
 0.0+0.0im  0.0+1.0im
```
While a hermitian conjugate of the same vector returns a `RowVector`
parametrized by the type `ConjArray`
```julia
julia> xc = [0.0, 1.0im]'
1×2 RowVector{Complex{Float64},ConjArray{Complex{Float64},1,Array{Complex{Float64},1}}}:
 0.0-0.0im  0.0-1.0im
```

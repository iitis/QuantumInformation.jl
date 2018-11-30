# Vectors and matrices in Julia

<!-- # Vectors and matrices in Julia

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

Values of variables `xt` and `xc` are views of the value of
variable `x`. When used the column and row vectors behave like bra and
kets, for example `xc*x` denotes the inner product of *bra*
`xc` and *ket* `x`, while `x*xc` denotes its outer
product resulting in a two-index array.

The linear algebra library in `Julia` provides standard operations on
matrices and vectors that are designed to take in to the account the types of
objects. -->

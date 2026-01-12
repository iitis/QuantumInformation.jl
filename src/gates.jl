export sx, sy, sz, qft, hadamard, grover

sx = ComplexF64[0 1; 1 0]
sy = ComplexF64[0 1im; -1im 0]
sz = ComplexF64[1 0; 0 -1]

export 𝕀, I
if VERSION < v"1.2"
  @deprecate 𝕀 LinearAlgebra.I
  (I::UniformScaling)(n::Integer) = Diagonal(fill(I.λ, n))
end
𝕀(::Type{T}, dim=2) where T<: Number = Matrix(one(T)*I(dim))
𝕀(dim=2) = 𝕀(ComplexF64, dim)


"""

- `d`: dimension of operator.

Prepares gate realized a [quantum Fourier transform](https://en.wikipedia.org/wiki/Quantum_Fourier_transform) of dimension `d`.
"""
qft(::Type{T}, d::Int) where T<:Number = [exp(T(2)*T(π)*1im*i*j/d) for i=0:d-1, j=0:d-1]/sqrt(real(T)(d))
qft(d::Int) = qft(ComplexF64, d)

"""

- `d`: dimension of operator.

Prepares [Grover operator](https://en.wikipedia.org/wiki/Grover%27s_algorithm) of dimension `d`.
"""
grover(::Type{T}, dim::Int) where T<:Number = ones(T,dim,dim)*2/dim - I
grover(dim::Int) = grover(ComplexF64, dim)

"""

- `d`: dimension of operator.

Prepares [Hadamard operator](https://en.wikipedia.org/wiki/Hadamard_transform) of dimension `d`.
"""
function hadamard(::Type{T}, dim::Int) where T<:Number
  if !ispow2(dim)
    throw(ArgumentError("Hadamard dim has to be power of 2"))
  end

  d=floor(log2(dim))
  H=one(T)/sqrt(one(T)*2)*[1 1;1 -1]
  mtx = 1
  for i=1:d
    mtx = mtx ⊗ H
  end
  return mtx
end
hadamard(dim::Int) = hadamard(ComplexF64, dim)



qft(d::Int) = [exp(2π*1im*i*j/d) for i=0:d-1, j=0:d-1]/sqrt(d)

grover(dim::Int) = ones(ComplexF64,dim,dim)*2/dim-diagm(ones(ComplexF64,dim))

function hadamard(dim::Int)
  if floor(log2(dim))!=log2(dim)
    throw(ArgumentError("Hadamard dim has to be power of 2"))
  end

  d=floor(log2(dim))
  H=1/sqrt(2)*[1 1;1 -1]
  mtx = 1
  for i=1:d
    mtx = mtx ⊗ H
  end
  return mtx
end

sx = ComplexF64[0 1; 1 0]
sy = ComplexF64[0 1im; -1im 0]
sz = ComplexF64[1 0; 0 -1]

function qft(dim::Int)
  mtx=zeros(ComplexF64,dim,dim)
  twopii = 2*pi*1im
  for i=0:dim-1
    for j=0:dim-1
      mtx[i+1,j+1]=exp(twopii*i*j/dim)
    end
  end
  mtx=mtx/sqrt(dim)
  return mtx
end

grover(dim::Int) = ones(ComplexF64,dim,dim)*2/dim-diagm(ones(ComplexF64,dim))

function hadamard(dim::Int)
  if(floor(log2(dim))!=log2(dim))
    error("hadamard: dim has to be power of 2")
  end
  d=floor(log2(dim))
  H=1/sqrt(2)*[1 1;1 -1]
  mtx=reduce(kron, [H for i=1:d])
  return mtx
end

sx = ComplexF64[0 1; 1 0]
sy = ComplexF64[0 1im; -1im 0]
sz = ComplexF64[1 0; 0 -1]

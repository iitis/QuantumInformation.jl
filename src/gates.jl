function qft(dim::Int)
  mtx=zeros(Complex128,dim,dim)
  twopii = 2*pi*1im
  for i=0:dim-1
    for j=0:dim-1
      mtx[i+1,j+1]=exp(twopii*i*j/dim)
    end
  end
  mtx=mtx/sqrt(dim)
  return mtx
end

grover(dim::Int) = ones(Complex128,dim,dim)*2/dim-diag(Complex128,ones(dim,1))

function hadamard(dim::Int)
  if(floor(log2(dim))!=log2(dim))
    error("hadamard: dim has to be power of 2")
  end
  d=floor(log2(dim))
  H=1/sqrt(2)*[1 1;1 -1]
  mtx=reduce(kron, [H for i=1:d])
  return mtx
end
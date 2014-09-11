ket(val, dim) = (k=zeros(Complex128,dim);k[val+1]=1.0;k)
bra(val, dim) = ket(val,dim)'
ketbra(valk, valb, dim) = (kb=zeros(Complex128,dim,dim);kb[valk+1,valb+1]=1.0;kb)
proj(ket) = ket*ket'
superoperator(kraus_list) = sum(k -> kron(k,k'), kraus_list)

random_ginibre_matrix(m,n) = G=randn(n,m)+im*randn(n,m)
function random_mixed_state_hs(d)
  A=random_ginibre_matrix(d,d)
  A=A*A'
  A=A/trace(A)
  return A
end

function random_ket(d)
  c=randn(d,1)+i*randn(d,1)
  c=c/norm(c)
  return c
end

function random_pure_state(d)
    return proj(random_ket(d))
end

base_matrices(dim) = @task for i=0:dim-1
  for j=0:dim-1
    produce(ketbra(i,j,dim))
  end
end

function channel_to_superoperator(channel,dim)
    Eijs=base_matrices(dim)
    A=[res(channel(e)) for e in Eijs]
    return A
end

res(X) = vec(permutedims(X,[2 1]))

function unres(X)
  s=sqrt(size(X,1))
  Xu=permutedims(reshape(X,s,s),invperm([2,1]))
end

applykraus(kraus_list,stin) = sum(k-> k*stin*k', kraus_list)

function ptrace(rho,idims,isystems)
    # convert notation to column-major form
    dims=reverse(idims)
    systems=length(idims)-isystems+1

    if size(rho,1)!=size(rho,2)
        error("Non square matrix passed to ptrace")
    end
    if prod(dims)!=size(rho,1)
        error("Product of dimensions do not match shape of matrix.")
    end
    if ! ((maximum(systems) <= length(dims) |  (minimum(systems) > length(dims))))
        error("System index out of range")
    end
    offset=length(dims)
    keep=setdiff(1:offset, systems)
    dispose=systems
    perm =[dispose,keep, dispose+offset,keep+offset]
    tensor=reshape(rho,[dims, dims]...)
    keepdim=prod([size(tensor,x) for x in keep])
    disposedim=prod([size(tensor,x) for x in dispose])
    tensor=permutedims(tensor,perm)

    tensor=reshape(tensor,disposedim,keepdim,disposedim,keepdim)
    ret = zeros(typeof(rho[1,1]),keepdim,keepdim)
    for i=1:keepdim
      for j=1:keepdim
#       ret[i,j]=sum(k->tensor[k,i,k,j], 1:disposedim)
        ret[i,j]=sum([tensor[k,i,k,j] for k in 1:disposedim])
      end
    end
    return ret
end

function number2mixedradix(n, bases)
  if n >= prod(bases)
    error("number to big to transform")
  end

  digits = Array(Int64,length(bases))
  i=1
  for base in reverse(bases)
      digits[i] = mod(n,base)
      n = div(n,base)
      i+=1
  end
  return digits
end

function mixedradix2number(digits, bases)
  if length(digits)>length(bases)
    error("more digits than radices")
  end
  res = 0
  for i=1:length(digits)
    res = res * bases[i] + digits[i]
  end
  return res
end

function qft(dim)
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

function grover(dim)
  return ones(Complex128,dim,dim)*2/dim-diag(Complex128,ones(dim,1))
end

function hadamard(dim)
  if(floor(log2(dim))!=log2(dim))
    error("hadamard: dim has to be power of 2")
  end
  d=floor(log2(dim))
  H=1/sqrt(2)*[1 1;1 -1]
  mtx=reduce(kron, [H for i=1:d])
  return mtx
end

function reshuffle(rho::Matrix)
  """
    Performs reshuffling of indices of a matrix.
    Given multiindexed matrix M_{(m,\mu),(n,\nu)} it returns
    matrix M_{(m,n),(\mu,\nu)}.
  """
  (r,c) = size(rho)
  sqrtr = int(sqrt(r))
  sqrtc = int(sqrt(c))
  dimrows = [sqrtr, sqrtr]
  dimcolumns = [sqrtc, sqrtc]
  tensor=reshape(rho, [dimrows, dimcolumns]...)
  perm=[1:4]
  (perm[2],perm[3])=(perm[3],perm[2])
  tensor=permutedims(tensor, perm)
  (r1,r2,c1,c2)=size(tensor)
  return reshape(tensor, (r1*r2,c1*c2)...)
end

include("utils.jl")

ket(val, dim) = (k=zeros(Complex128,dim);k[val+1]=1.0;k)

bra(val, dim) = ket(val,dim)'

ketbra(valk, valb, dim) = (kb=zeros(Complex128,dim,dim);kb[valk+1,valb+1]=1.0;kb)

proj(ket) = ket*ket'

base_matrices(dim) = [ketbra(i,j,dim)for i=0:dim-1, j=0:dim-1]

res(X) = vec(permutedims(X,[2 1]))

function unres(X)
  s=int(sqrt(size(X,1)))
  Xu=permutedims(reshape(X,s,s),invperm([2,1]))
end

function kraus_to_superoperator(kraus_list)
  return sum((k) -> kron(k,k'), kraus_list)
end

function channel_to_superoperator(channel,dim)
    Eijs=base_matrices(dim)
    A=[res(channel(e)) for e in Eijs]
    return A
end

apply_kraus(kraus_list,stin) = sum(k-> k*stin*k', kraus_list)

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
  perm = [4, 2, 3, 1]
  tensor=permutedims(tensor, perm)
  (r1,r2,c1,c2)=size(tensor)
  return reshape(tensor, (r1*r2,c1*c2)...)
end

trace_distance(ρ, σ) = sum(abs(eigvals(Hermitian(ρ - σ))))

function fidelity_sqrt(ρ, σ)
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    error("Non square matrix detected")
  end
  vals = real(eigvals(ρ * σ))
  return sum(sqrt(real(vals[vals.>0])))
end

function fidelity(ρ, σ)
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    error("Non square matrix detected")
  end
  return fidelity_sqrt(ρ, σ)^2
end

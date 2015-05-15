using Devectorize
include("utils.jl")

function ket{T<:Union(Float64, Complex128)}(M::Type{T}, val::Int64, dim::Int64)
    ϕ=zeros(M, dim)
    ϕ[val+1]=1.0
    ϕ
end

ket(val::Int64, dim::Int64) = ket(Complex128, val, dim)
bra{T<:Union(Float64, Complex128)}(M::Type{T}, val::Int64, dim::Int64) = ket(M, val, dim)'
bra(val::Int64, dim::Int64) = bra(Complex128, val, dim)

function ketbra{T<:Union(Float64, Complex128)}(M::Type{T}, valk::Int64, valb::Int64, dim::Int64)
    ϕψ=zeros(M, dim, dim)
    ϕψ[valk+1,valb+1]=1.0
    ϕψ
end

ketbra(valk::Int64, valb::Int64, dim::Int64) = ketbra(Complex128, valk, valb, dim)

proj{T<:Union(Float64, Complex128)}(ket::Vector{T}) = ket * ket'

base_matrices(dim) = [ketbra(j, i, dim) for i=0:dim-1, j=0:dim-1] #TODO: should be a generator, but sometimes causes errors in applications

res{T<:Union(Float64, Complex128)}(ρ::Matrix{T}) = vec(permutedims(ρ, [2 1]))

#TODO: verify whether invperm is squivalent to transpose in this case
#TODO: user array views for reshaping?
function unres{T<:Union(Float64, Complex128)}(ϕ::Vector{T})
    s=int(sqrt(size(ϕ, 1)))
    permutedims(reshape(ϕ, s, s),invperm([2,1]))
end

unres{T<:Union(Float64, Complex128)}(ϕ::Vector{T}, m::Int64, n::Int64) = permutedims(reshape(ϕ, n, m), invperm([2,1]))

kraus_to_superoperator{T<:Union(Float64, Complex128)}(kraus_list::Vector{Matrix{T}}) = sum((k) -> kron(k, k'), kraus_list)

function channel_to_superoperator(channel::Function, dim::Int64)
    #TODO: create a matrix of zeros first, then only set columns
    Eijs=base_matrices(dim)
    hcat([res(channel(e)) for e in Eijs]...)
end

apply_kraus{T<:Union(Float64, Complex128)}(kraus_list::Vector{Matrix{T}}, ρ::Matrix{T}) = sum(k-> k*ρ*k', kraus_list)

function ptrace{T<:Union(Float64, Complex128)}(ρ::Matrix{T}, idims::Vector, isystems::Vector)
    # TODO: convert notation to column-major form
    dims=reverse(idims)
    systems=length(idims)-isystems+1

    if size(ρ,1)!=size(ρ,2)
        error("Non square matrix passed to ptrace")
    end
    if prod(dims)!=size(ρ,1)
        error("Product of dimensions do not match shape of matrix.")
    end
    if ! ((maximum(systems) <= length(dims) |  (minimum(systems) > length(dims))))
        error("System index out of range")
    end
    offset=length(dims)
    keep=setdiff(1:offset, systems)
    dispose=systems
    perm =[dispose,keep, dispose+offset,keep+offset]
    tensor=reshape(ρ,[dims, dims]...)
    keepdim=prod([size(tensor,x) for x in keep])
    disposedim=prod([size(tensor,x) for x in dispose])
    tensor=permutedims(tensor,perm)

    tensor=reshape(tensor,disposedim,keepdim,disposedim,keepdim)
    ret = zeros(typeof(ρ[1,1]),keepdim,keepdim)
    for i=1:keepdim
      for j=1:keepdim
        ret[i,j]=sum([tensor[k,i,k,j] for k in 1:disposedim])
      end
    end
    return ret
end
#TODO: write tests for thuis
#TODO: allow for more than bipartite systems???
function ptrace{T<:Union(Float64, Complex128)}(ϕ::Vector{T}, idims::Vector, isystem::Int64)
    A = unres(ϕ, idims...)
    if isystem == 1
        return A'*A
    elseif isystem == 2
        return A*A'
    end
end

function number2mixedradix(n::Int64, bases::Vector{Int64})
    if n >= prod(bases)
        error("number to big to transform")
    end

    digits = Array(Int64, length(bases))
    for (i, base) in enumerate(reverse(bases))
        n, digits[i] = divrem(n, base)
    end
    digits
end

function mixedradix2number(digits::Vector{Int64}, bases::Vector{Int64})
    if length(digits)>length(bases)
        error("more digits than radices")
    end
    
    res = 0
    digitsreversed = reverse(digits)
    for i=1:length(digits)
        res = res * bases[i] + digitsreversed[i]
    end
    res
end

function reshuffle{T<:Union(Float64, Complex128)}(ρ::Matrix{T})
  """
    Performs reshuffling of indices of a matrix.
    Given multiindexed matrix M_{(m,\mu),(n,\nu)} it returns
    matrix M_{(m,n),(\mu,\nu)}.
  """
  (r, c) = size(ρ)
  sqrtr = int(sqrt(r))
  sqrtc = int(sqrt(c))
  dimrows = [sqrtr, sqrtr]
  dimcolumns = [sqrtc, sqrtc]
  tensor=reshape(ρ, [dimrows, dimcolumns]...)
  perm = [4, 2, 3, 1]
  tensor=permutedims(tensor, perm)
  (r1,r2,c1,c2)=size(tensor)
  return reshape(tensor, (r1*r2,c1*c2)...)
end

trace_distance{T<:Union(Float64, Complex128)}(ρ::Matrix{T}, σ::Matrix{T}) = sum(abs(eigvals(Hermitian(ρ - σ))))

function fidelity_sqrt{T<:Union(Float64, Complex128)}(ρ::Matrix{T}, σ::Matrix{T})
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    error("Non square matrix detected")
  end
  λ = real(eigvals(ρ * σ))
  @devec r = sum(sqrt(λ[λ.>0]))
end

function fidelity{T<:Union(Float64, Complex128)}(ρ::Matrix{T}, σ::Matrix{T})
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    error("Non square matrix detected")
  end
  return fidelity_sqrt(ρ, σ)^2
end

fidelity{T<:Union(Float64, Complex128)}(ϕ::Vector{T}, ψ::Vector{T}) = abs2(dot(ϕ, ψ))
fidelity{T<:Union(Float64, Complex128)}(ϕ::Vector{T}, ρ::Matrix{T}) = ϕ' * ρ * ϕ
fidelity{T<:Union(Float64, Complex128)}(ρ::Matrix{T}, ϕ::Vector{T}) = fidelity(ϕ, ρ)

function shannon_entropy{T<:Real}(p::Vector{T})
    @devec r = -sum(p .* log(p))
end

shannon_entropy{T<:Real}(x::T) = x > 0 ? -x * log(x) - (1 - x) * log(1 - x) : error("Negative number passed to shannon_entropy")

function entropy{T<:Union(Float64, Complex128)}(ρ::Hermitian{T})
    λ = eigvals(ρ)
    λ = λ[λ .> 0]
    @devec r = -sum(λ .* log(λ))
end

entropy{T<:Union(Float64, Complex128)}(H::Matrix{T}) = ishermitian(H) ? entropy(Hermitian(H)) : error("Non-hermitian matrix passed to entropy")
entropy{T<:Union(Float64, Complex128)}(ϕ::Vector{T}) = 0.

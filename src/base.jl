function ket{T<:Union{Float64, Complex128}}(::Type{T}, val::Int64, dim::Int64)
    ϕ=zeros(T, dim)
    ϕ[val+1]=1.0
    ϕ
end

ket(val::Int64, dim::Int64) = ket(Complex128, val, dim)
bra{T<:Union{Float64, Complex128}}(::Type{T}, val::Int64, dim::Int64) = ket(T, val, dim)'
bra(val::Int64, dim::Int64) = bra(Complex128, val, dim)

function ketbra{T<:Union{Float64, Complex128}}(::Type{T}, valk::Int64, valb::Int64, dim::Int64)
    ϕψ=zeros(T, dim, dim)
    ϕψ[valk+1,valb+1]=1.0
    ϕψ
end

ketbra(valk::Int64, valb::Int64, dim::Int64) = ketbra(Complex128, valk, valb, dim)

proj{T<:Union{Float64, Complex128}}(ket::Vector{T}) = ket * ket'

function base_matrices(dim)
    function _it()
        for i=0:dim-1, j=0:dim-1
            produce(ketbra(j, i, dim))
        end
    end
    Task(_it)
end

res{T<:Union{Float64, Complex128}}(ρ::Matrix{T}) = vec(transpose(ρ))

function unres{T<:Union{Float64, Complex128}}(ϕ::Vector{T})
    s=round(Int, sqrt(size(ϕ, 1)), RoundDown)
    transpose(reshape(ϕ, s, s))
end

unres{T<:Union{Float64, Complex128}}(ϕ::Vector{T}, m::Int64, n::Int64) = transpose(reshape(ϕ, n, m))

kraus_to_superoperator{T<:Union{Float64, Complex128}}(kraus_list::Vector{Matrix{T}}) = sum((k) -> kron(k, k'), kraus_list)

function channel_to_superoperator(channel::Function, dim::Int64)
    M = zeros(Complex128, dim*dim, dim*dim)
    for (i, e) in enumerate(base_matrices(dim))
        M[:, i] = res(channel(e))
    end
    M
end

apply_kraus{T<:Union{Float64, Complex128}}(kraus_list::Vector{Matrix{T}}, ρ::Matrix{T}) = sum(k-> k*ρ*k', kraus_list)

function ptrace{T<:Union{Float64, Complex128}}(ρ::Matrix{T}, idims::Vector, isystems::Vector)
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
    offset = length(dims)
    keep = setdiff(1:offset, systems)
    dispose = systems
    perm  = [dispose; keep; dispose+offset; keep+offset]
    tensor = reshape(ρ, [dims; dims]...)
    keepdim = prod([size(tensor,x) for x in keep])
    disposedim = prod([size(tensor,x) for x in dispose])
    tensor = permutedims(tensor,perm)

    tensor=reshape(tensor, disposedim, keepdim, disposedim, keepdim)
    ret = zeros(typeof(ρ[1,1]),keepdim,keepdim)
    for i=1:keepdim
      for j=1:keepdim
        ret[i,j]=sum([tensor[k,i,k,j] for k in 1:disposedim])
      end
    end
    return ret
end
#TODO: allow for more than bipartite systems???
function ptrace{T<:Union{Float64, Complex128}}(ϕ::Vector{T}, idims::Vector, isystem::Int64)
    A = unres(ϕ, idims...)
    if isystem == 1
        return A'*A
    elseif isystem == 2
        return A*A'
    end
end

function ptranspose{T<:Union{Float64, Complex128}}(ρ::Matrix{T}, idims::Vector, isystems::Vector)
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

    offset = length(dims)
    tensor = reshape(ρ, [dims; dims]...)
    perm = collect(1:(2offset))
    for s in systems
        idx1 = find(x->x==s, perm)[1]
        idx2 = find(x->x==(s + offset), perm)[1]
        perm[idx1], perm[idx2] = perm[idx2], perm[idx1]
    end
    tensor = permutedims(tensor, invperm(perm))
    reshape(tensor, size(ρ))
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
# FIX THESE
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

function reshuffle{T<:Union{Float64, Complex128}}(ρ::Matrix{T})
  """
    Performs reshuffling of indices of a matrix.
    Given multiindexed matrix M_{(m,\mu),(n,\nu)} it returns
    matrix M_{(m,n),(\mu,\nu)}.
  """
  (r, c) = size(ρ)
  sqrtr = round(Int, sqrt(r), RoundDown)
  sqrtc = round(Int, sqrt(c), RoundDown)
  tensor = reshape(ρ, sqrtr, sqrtr, sqrtc, sqrtc)
  perm = [4, 2, 3, 1]
  tensor = permutedims(tensor, perm)
  (r1, r2, c1, c2) = size(tensor)
  return reshape(tensor, r1*r2, c1*c2)
end

trace_distance{T<:Union{Float64, Complex128}}(ρ::Matrix{T}, σ::Matrix{T}) = sum(abs(eigvals(Hermitian(ρ - σ))))

function fidelity_sqrt{T<:Union{Float64, Complex128}}(ρ::Matrix{T}, σ::Matrix{T})
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    error("Non square matrix detected")
  end
  λ = real(eigvals(ρ * σ))
  r = sum(sqrt.(λ[λ.>0]))
end

function fidelity{T<:Union{Float64, Complex128}}(ρ::Matrix{T}, σ::Matrix{T})
  if size(ρ, 1) != size(ρ, 2) || size(σ, 1) != size(σ, 2)
    error("Non square matrix detected")
  end
  return fidelity_sqrt(ρ, σ)^2
end

fidelity{T<:Union{Float64, Complex128}}(ϕ::Vector{T}, ψ::Vector{T}) = abs2(dot(ϕ, ψ))
fidelity{T<:Union{Float64, Complex128}}(ϕ::Vector{T}, ρ::Matrix{T}) = ϕ' * ρ * ϕ
fidelity{T<:Union{Float64, Complex128}}(ρ::Matrix{T}, ϕ::Vector{T}) = fidelity(ϕ, ρ)

shannon_entropy{T<:Real}(p::Vector{T}) = -sum(p .* log.(p))

shannon_entropy{T<:Real}(x::T) = x > 0 ? -x * log(x) - (1 - x) * log(1 - x) : error("Negative number passed to shannon_entropy")

function entropy{T<:Union{Float64, Complex128}}(ρ::Hermitian{T})
    λ = eigvals(ρ)
    λ = λ[λ .> 0]
    -sum(λ .* log(λ))
end

entropy{T<:Union{Float64, Complex128}}(H::Matrix{T}) = ishermitian(H) ? entropy(Hermitian(H)) : error("Non-hermitian matrix passed to entropy")
entropy{T<:Union{Float64, Complex128}}(ϕ::Vector{T}) = 0.

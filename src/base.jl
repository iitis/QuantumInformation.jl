function ket(::Type{Tv}, val::Int, dim::Int) where Tv<:AbstractVector{T} where T<:Number
    dim > 0 ? () : error("Vector dimension has to be nonnegative")
    val < dim ? () : error("Label have to be smaller than vector dimmension")
    ϕ = zeros(T, dim)
    ϕ[val+1] = one(T)
    ϕ
end

function ket(::Type{Tv}, val::Int, dim::Int) where Tv<:AbstractSparseVector{T} where T<:Number
    dim > 0 ? () : error("Vector dimension has to be nonnegative")
    val < dim ? () : error("Label have to be smaller than vector dimmension")
    ϕ = spzeros(T, dim)
    ϕ[val+1] = one(T)
    ϕ
end

function ket(val::Int, dim::Int; sparse=false)
    if sparse==true
        return ket(SparseVector{ComplexF64}, val, dim)
    else
        return ket(Vector{ComplexF64}, val, dim)
    end
end

bra(::Type{Tv}, val::Int, dim::Int) where Tv<:AbstractVector{T} where T<:Number = ket(Tv, val, dim)'

function bra(val::Int, dim::Int; sparse=false)
    if sparse==true
        return bra(SparseVector{ComplexF64}, val, dim)
    else
        return bra(Vector{ComplexF64}, val, dim)
    end
end

function ketbra(::Type{Tv}, valk::Int, valb::Int, dim::Int) where Tv<:AbstractMatrix{T} where T<:Number
    dim > 0 ? () : error("Vector dimension has to be nonnegative")
    valk < dim && valb < dim ? () : error("Ket and bra labels have to be smaller than operator dimmension")
    ϕψ = zeros(T, dim, dim)
    ϕψ[valk+1,valb+1] = one(T)
    ϕψ
end

function ketbra(::Type{Tv}, valk::Int, valb::Int, dim::Int) where Tv<:AbstractSparseMatrix{T} where T<:Number
    dim > 0 ? () : error("Vector dimension has to be nonnegative")
    valk < dim && valb < dim ? () : error("Ket and bra labels have to be smaller than operator dimmension")
    ϕψ = spzeros(T, dim, dim)
    ϕψ[valk+1,valb+1] = one(T)
    ϕψ
end

function ketbra(valk::Int, valb::Int, dim::Int; sparse=false)
    if sparse==true
        return ketbra(SparseMatrixCSC{ComplexF64}, valk, valb, dim)
    else
        return ketbra(Matrix{ComplexF64}, valk, valb, dim)
    end
end

proj(ket::AbstractVector{T}) where T<:Number = ket * ket'

# function base_matrices(dim)
#     function _it()
#         for i=0:dim-1, j=0:dim-1
#             produce(ketbra(j, i, dim))
#         end
#     end
#     Task(_it)
# end

base_matrices(dim) = Channel() do c
    dim > 0 ? () : error("Operator dimension has to be nonnegative")
    for i=0:dim-1, j=0:dim-1
        push!(c, ketbra(j, i, dim))
    end
end

res(ρ::AbstractMatrix{T}) where T<:Number = vec(transpose(ρ))

function unres(ϕ::AbstractVector{T}) where T<:Number
    dim = size(ϕ, 1)
    s = isqrt(dim)
    s*s == dim ? () : error("Vector has to have dimensions being a square of an integer number")
    transpose(reshape(ϕ, s, s))
end

unres(ϕ::AbstractVector{T}, m::Int, n::Int) where T<:Number = transpose(reshape(ϕ, n, m))

# TODO: allow different type of Kraus operators and the quantum state
function apply_kraus(kraus_list::Vector{T}, ρ::T) where {T<:AbstractMatrix{T1}} where {T1<:Number}
    # TODO: chceck if all Kraus operators are the same shape and fit the input state
    sum(k-> k*ρ*k', kraus_list)
end

function ptrace(ρ::AbstractMatrix{T}, idims::Vector, isystems::Vector) where T<:Number
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

# TODO: allow for more than bipartite systems???
function ptrace(ϕ::AbstractVector{T}, idims::Vector, isystem::Int) where T<:Number
    A = unres(ϕ, idims...)
    if isystem == 1
        return A'*A
    elseif isystem == 2
        return A*A'
    end
end

function ptranspose(ρ::AbstractMatrix{T}, idims::Vector, isystems::Vector) where T<:Number
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

function number2mixedradix(n::Int, bases::Vector{Int})
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
function mixedradix2number(digits::Vector{Int}, bases::Vector{Int})
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

"""
  Performs reshuffling of indices of a matrix.
  Given multiindexed matrix M_{(m,μ),(n,ν)} it returns
  matrix M_{(m,n),(μ,ν)}.
"""
function reshuffle(ρ::AbstractMatrix{T}) where T<:Number
  (r, c) = size(ρ)
  sqrtr = isqrt(r)
  sqrtc = isqrt(c)
  tensor = reshape(ρ, sqrtr, sqrtr, sqrtc, sqrtc)
  perm = [4, 2, 3, 1]
  tensor = permutedims(tensor, perm)
  (r1, r2, c1, c2) = size(tensor)
  return reshape(tensor, r1*r2, c1*c2)
end

function max_mix(dim):
    return (1.0 / dim * eye(dim))
end

function max_entangled(dim):
    sqrtdim = isqrt(dim)
    1 / sqrt(sqrtdim) * sum(kron(ket(i, sqrtdim), ket(i, sqrtdim)) for i in 1:sqrtdim)
end

"""
http://en.wikipedia.org/wiki/Werner_state
"""
function werner_state(alpha, dim):
    return alpha * proj(max_entangled(dim)) + (1 - alpha) * max_mix(dim)
end

#=
TODO: port to julia
def base_hermitian_matrices(dim):
    """
    Generator. Returns elementary hermitian matrices of dimension dim x dim.
    """
    for (a, b) in product(xrange(dim), repeat=2):
        if a > b:
            yield 1 / np.sqrt(2) * np.matrix(1j * ketbra(a, b, dim) - 1j * ketbra(b, a, dim))
        elif a < b:
            yield 1 / np.sqrt(2) * np.matrix(ketbra(a, b, dim) + ketbra(b, a, dim))
        else:
            yield np.matrix(ketbra(a, b, dim))


def permute_systems(rho, dims, systemperm):
    rho = np.asarray(rho)
    dims = list(dims)
    systemperm = list(systemperm)
    if rho.shape[0] != rho.shape[1]:
        raise Exception("Non square matrix passed to ptrace")
    if np.prod(dims) != rho.shape[0]:
        raise Exception("Product of dimensions do not match shape of matrix.")
    if not ((max(systemperm) <= len(dims) or (min(systemperm) > len(dims)))):
        raise Exception("System index out of range")
    offset = len(dims)
    perm1 = systemperm
    perm2 = map(lambda x: x + offset, perm1)
    perm = perm1 + perm2
    tensor = np.array(rho).reshape(2 * dims)
    tensor = tensor.transpose(perm)
    return np.asmatrix(tensor.reshape(rho.shape))
=#

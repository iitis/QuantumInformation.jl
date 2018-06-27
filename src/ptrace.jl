
function ptrace(ρ::AbstractMatrix{<:Number}, idims::Vector{Int}, isystems::Vector{Int})
    # TODO: convert notation to column-major form
    dims=reverse(idims)
    systems=length(idims)-isystems+1

    if size(ρ,1)!=size(ρ,2)
        throw(ArgumentError("Non square matrix passed to ptrace"))
    end
    if prod(dims)!=size(ρ,1)
        throw(ArgumentError("Product of dimensions do not match shape of matrix."))
    end
    if maximum(systems) > length(dims) || minimum(systems) < 1
        throw(ArgumentError("System index out of range"))
    end
    offset = length(dims)
    keep = setdiff(1:offset, systems)
    dispose = systems
    perm  = [dispose; keep; dispose+offset; keep+offset]
    tensor = reshape(ρ, [dims; dims]...)
    keepdim = prod([size(tensor, x) for x in keep])
    disposedim = prod([size(tensor, x) for x in dispose])
    tensor = permutedims(tensor, perm)

    tensor=reshape(tensor, disposedim, keepdim, disposedim, keepdim)
    ret = zeros(typeof(ρ[1,1]), keepdim, keepdim)
    for i=1:keepdim
      for j=1:keepdim
        ret[i,j] = sum([tensor[k, i, k, j] for k in 1:disposedim]) # TODO: change to iterator
      end
    end
    return ret
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `idims`: dimensins of subsystems.
- `sys`: traced subsystem.

Return [partial trace](https://en.wikipedia.org/wiki/Partial_trace) of matrix `ρ` over the subsystems determined by `isystems`.
"""
ptrace(ρ::AbstractMatrix{<:Number}, idims::Vector{Int}, sys::Int) = ptrace(ρ, idims, [sys])

# TODO: allow for more than bipartite systems???
function ptrace(ϕ::AbstractVector{<:Number}, idims::Vector{Int}, sys::Int)
    _, cols = idims
    m = unres(ϕ, cols)
    length(idims) == 2 ? () : throw(ArgumentError("idims has to be of length 2"))
    if sys == 1
        return m'*m
    elseif sys == 2
        return m*m'
    else
        throw(ArgumentError("sys must be 1 or 2"))
    end
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `idims`: dimensins of subsystems.
- `isystems`: traced subsystems.

Return [partial trace](https://en.wikipedia.org/wiki/Partial_trace) of matrix `ρ` over the subsystems determined by `isystems`.
"""
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

ptrace(ρ::AbstractMatrix{<:Number}, idims::Vector{Int}, sys::Int) = ptrace(ρ, idims, [sys])

function ptrace(ρ::AbstractSparseMatrix{T}, idims::Vector{Int}, sys::Int) where T<:Number

    if size(ρ,1)!=size(ρ,2)
        throw(ArgumentError("Non square matrix passed to ptrace"))
    end
    if prod(idims)!=size(ρ,1)
        throw(ArgumentError("Product of dimensions do not match shape of matrix."))
    end
    if sys > 2 || sys < 1
        throw(ArgumentError("sys must be either 1 or 2, not $sys"))
    end
    if length(idims) != 2
        throw(ArgumentError("Only bipartite systems supported"))
    end

    d1, d2 = idims
    if sys == 1
        return sum([ρ[i*d2+1:(i+1)*d2, i*d2+1:(i+1)*d2] for i=0:d1-1])
    elseif sys == 2
        I, J, V = Int[], Int[], T[]
        for i=0:d1-1, j=0:d1-1
            v = trace(ρ[i*d2+1:(i+1)*d2, j*d2+1:(j+1)*d2])
            if isapprox(v, 0)
                continue
            end
            push!(I, i+1)
            push!(J, j+1)
            push!(V, v)
        end
        return sparse(I, J, V, d1, d1)
    end
end

ptrace(ρ::AbstractSparseMatrix{<:Number}, idims::Vector{Int}, sys::Vector{Int}) = ptrace(ρ, idims, sys[1])

# TODO: allow for more than bipartite systems???
function ptrace(ϕ::AbstractVector{<:Number}, idims::Vector{Int}, sys::Int)
    _, cols = idims
    A = unres(ϕ, cols)
    if sys == 1
        return A'*A
    elseif sys == 2
        return A*A'
    else
        throw(ArgumentError("sys must be 1 or 2"))
    end

end

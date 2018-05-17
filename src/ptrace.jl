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

ptrace(ρ::AbstractMatrix{T}, idims::Vector, sys::Int) where T<:Number = ptrace(ρ, idims, [sys])

function ptrace(ρ::AbstractSparseMatrix{T}, idims::Vector, sys::Int) where T<:Number

    if size(ρ,1)!=size(ρ,2)
        error("Non square matrix passed to ptrace")
    end
    if prod(idims)!=size(ρ,1)
        error("Product of dimensions do not match shape of matrix.")
    end
    if sys > 2 || sys < 1
        error("sys mus be either 1 or 2, not $sys")
    end
    if length(idims) != 2
        error("Only bipartite systems supported")
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

# TODO: allow for more than bipartite systems???
function ptrace(ϕ::AbstractVector{T}, idims::Vector, isystem::Int) where T<:Number
    A = unres(ϕ, idims...)
    if isystem == 1
        return A'*A
    elseif isystem == 2
        return A*A'
    end
end

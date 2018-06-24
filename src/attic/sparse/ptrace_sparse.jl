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

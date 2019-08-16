export ptrace

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `idims`: dimensins of subsystems.
- `isystems`: traced subsystems.

Return [partial trace](https://en.wikipedia.org/wiki/Partial_trace) of matrix `ρ` over the subsystems determined by `isystems`.
"""
function ptrace(ρ::AbstractMatrix{<:Number}, idims::Vector{Int}, isystems::Vector{Int})
    dims = reverse(idims)
    systems = length(idims) .- isystems .+ 1

    if size(ρ,1) != size(ρ,2)
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

    traceidx = [1:offset; 1:offset]
    traceidx[keep] .+= offset

    tensor = reshape(ρ, [dims; dims]...)
    keepdim = prod([size(tensor, x) for x in keep])
    return reshape(tensortrace(tensor, Tuple(traceidx)), keepdim, keepdim)
end

"""
$(SIGNATURES)
- `ρ`: quantum state.
- `idims`: dimensins of subsystems.
- `sys`: traced subsystem.
"""
ptrace(ρ::AbstractMatrix{<:Number}, idims::Vector{Int}, sys::Int) = ptrace(ρ, idims, [sys])

"""
$(SIGNATURES)
- `ψ`: quantum state pure state (ket).
- `idims`: dimensins of subsystems - only bipartite states accepted.
- `sys`: traced subsystem.
"""
function ptrace(ψ::AbstractVector{<:Number}, idims::Vector{Int}, sys::Int)
    # TODO : Allow mutlipartite systems
    _, cols = idims
    m = unres(ψ, cols)
    length(idims) == 2 ? () : throw(ArgumentError("idims has to be of length 2"))
    if sys == 1
        return m'*m
    elseif sys == 2
        return m*m'
    else
        throw(ArgumentError("sys must be 1 or 2"))
    end
end

export ptrace
using SparseArrays

"""
  - `ρ`: quantum state.
  - `idims`: dimensins of subsystems.
  - `isystems`: traced subsystems.

Return [partial trace](https://en.wikipedia.org/wiki/Partial_trace) of matrix `ρ` over the subsystems determined by `isystems`.
"""
function ptrace(ρ::AbstractMatrix, idims::Vector{Int}, isystems::Vector{Int})
    dims = reverse(idims)
    systems = length(idims) .- isystems .+ 1

    if size(ρ, 1) != size(ρ, 2)
        throw(ArgumentError("Non square matrix passed to ptrace"))
    end
    if prod(dims)!=size(ρ, 1)
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

function ptrace(ρ::AbstractSparseMatrix, idims::Vector{Int}, isystems::Vector{Int})
    size(ρ, 1) != size(ρ, 2) && throw(ArgumentError("Non-square matrix"))

    dims = reverse(idims)
    systems = length(idims) .- isystems .+ 1

    prod(dims) != size(ρ, 1) && throw(ArgumentError("Dimensions mismatch"))

    n_sys = length(dims)
    # Identify systems to keep and trace
    keep_indices = setdiff(1:n_sys, systems)

    # Calculate strides for index decomposition
    strides = [1]
    for d in dims[1:(end - 1)]
        push!(strides, strides[end] * d)
    end

    # Pre-calculate dimensions for result
    keep_dims = dims[keep_indices]
    keep_dim_total = prod(keep_dims)

    I, J, V = findnz(ρ)

    # Dictionary to aggregate values for result
    res_dict = Dict{Tuple{Int, Int}, eltype(V)}()

    # Output strides
    out_strides = [1]
    for d in keep_dims[1:(end - 1)]
        push!(out_strides, out_strides[end] * d)
    end

    for k in 1:length(V)
        val = V[k]
        r = I[k] - 1
        c = J[k] - 1

        # Decompose indices
        r_digits = zeros(Int, n_sys)
        c_digits = zeros(Int, n_sys)

        curr_r = r
        curr_c = c

        match = true

        for i in 1:n_sys
            r_digits[i] = curr_r % dims[i]
            curr_r ÷= dims[i]

            c_digits[i] = curr_c % dims[i]
            curr_c ÷= dims[i]
        end

        # Verify trace systems match
        for sys in systems
            if r_digits[sys] != c_digits[sys]
                match = false
                break
            end
        end

        if match
            # Calculate new indices
            new_r = 0
            new_c = 0

            for (idx, sys) in enumerate(keep_indices)
                new_r += r_digits[sys] * out_strides[idx]
                new_c += c_digits[sys] * out_strides[idx]
            end

            key = (new_r + 1, new_c + 1)
            res_dict[key] = get(res_dict, key, zero(eltype(V))) + val
        end
    end

    II = Int[]
    JJ = Int[]
    VV = eltype(V)[]

    for ((r, c), v) in res_dict
        if v != 0
            push!(II, r)
            push!(JJ, c)
            push!(VV, v)
        end
    end

    return sparse(II, JJ, VV, keep_dim_total, keep_dim_total)
end

"""
  - `ρ`: quantum state.
  - `idims`: dimensins of subsystems.
  - `sys`: traced subsystem.
"""
ptrace(ρ::AbstractMatrix, idims::Vector{Int}, sys::Int) = ptrace(ρ, idims, [sys])

"""
  - `ψ`: quantum state pure state (ket).
  - `idims`: dimensins of subsystems - only bipartite states accepted.
  - `sys`: traced subsystem.
"""
function ptrace(ψ::AbstractVector, idims::Vector{Int}, sys::Int)
    # TODO : Allow mutlipartite systems
    length(idims) == 2 ? () : throw(ArgumentError("idims has to be of length 2"))
    _, cols = idims
    m = unres(ψ, cols)
    if sys == 1
        return transpose(m) * conj.(m)
    elseif sys == 2
        return m * m'
    else
        throw(ArgumentError("sys must be 1 or 2"))
    end
end

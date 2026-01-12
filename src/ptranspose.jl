export ptranspose
using SparseArrays
"""

- `ρ`: quantum state.
- `idims`: dimensins of subsystems.
- `isystems`: transposed subsystems.

Return [partial transposition](http://en.wikipedia.org/wiki/Peres-Horodecki_criterion) of matrix `ρ` over the subsystems determined by `isystems`.
"""
function ptranspose(ρ::AbstractMatrix, idims::Vector{Int}, isystems::Vector{Int})
    dims = reverse(idims)
    systems = length(idims) .- isystems .+ 1

    if size(ρ,1)!=size(ρ,2)
        throw(ArgumentError("Non square matrix passed to ptrace"))
    end
    if prod(dims)!=size(ρ,1)
        throw(ArgumentError("Product of dimensions do not match shape of matrix."))
    end
    if maximum(systems) > length(dims) ||  minimum(systems) < 1
        throw(ArgumentError("System index out of range"))
    end

    offset = length(dims)
    tensor = reshape(ρ, [dims; dims]...)
    perm = collect(1:(2offset))
    for s in systems
        idx1 = findfirst(x->x==s, perm)
        idx2 = findfirst(x->x==(s + offset), perm)
        perm[idx1], perm[idx2] = perm[idx2], perm[idx1]
    end
    tensor = permutedims(tensor, invperm(perm))
    reshape(tensor, size(ρ))
end

function ptranspose(ρ::AbstractSparseMatrix, idims::Vector{Int}, isystems::Vector{Int})
    size(ρ, 1) != size(ρ, 2) && throw(ArgumentError("Non-square matrix"))
    
    dims = reverse(idims)
    systems = length(idims) .- isystems .+ 1
    
    prod(dims) != size(ρ, 1) && throw(ArgumentError("Dimensions mismatch"))

    n_sys = length(dims)
    
    # Calculate strides
    strides = [1]
    for d in dims[1:end-1]
        push!(strides, strides[end] * d)
    end
    
    I, J, V = findnz(ρ)
    
    new_I = copy(I)
    new_J = copy(J)
    
    for k in 1:length(V)
        r = I[k] - 1
        c = J[k] - 1
        
        r_digits = zeros(Int, n_sys)
        c_digits = zeros(Int, n_sys)
        
        curr_r = r
        curr_c = c
        
        for i in 1:n_sys
            r_digits[i] = curr_r % dims[i]
            curr_r ÷= dims[i]
            
            c_digits[i] = curr_c % dims[i]
            curr_c ÷= dims[i]
        end
        
        # Swap indices for transposed systems
        for sys in systems
             r_digits[sys], c_digits[sys] = c_digits[sys], r_digits[sys]
        end
        
        # Reconstruct indices
        new_r = 0
        new_c = 0
        for i in 1:n_sys
            new_r += r_digits[i] * strides[i]
            new_c += c_digits[i] * strides[i]
        end
        
        new_I[k] = new_r + 1
        new_J[k] = new_c + 1
    end
    
    sparse(new_I, new_J, V, size(ρ, 1), size(ρ, 2))
end

"""

- `ρ`: quantum state.
- `idims`: dimensins of subsystems.
- `sys`: transposed subsystem.
"""
ptranspose(ρ::AbstractMatrix, idims::Vector{Int}, sys::Int) = ptranspose(ρ, idims, [sys])

function _ptranspose(ρ::AbstractMatrix{<:Number}, idims::Vector{Int}, isystems::Vector{Int})
    ns = length(idims)

    ex1 = Expr(:ref, :x)
    ex2 = Expr(:ref, ρ)

    I = Expr(:tuple, [gensym() for _=1:ns]...)
    J = Expr(:tuple, [gensym() for _=1:ns]...)

    K = copy(I)
    L = copy(J)

    r = Expr(:tuple)
    for (k, (i, j)) in enumerate(zip(K.args, L.args))
        push!(r.args, :($i in 1:$(idims[k])), :($j in 1:$(idims[k])))
    end
    for s in isystems
        K.args[s], L.args[s] = L.args[s], K.args[s]
    end
    push!(ex1.args, I, J)
    push!(ex2.args, L, K)

    ex = Expr(:(:=), ex1, ex2)
    ex, r
    @eval @cast $ex $r
end
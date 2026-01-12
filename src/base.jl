export ket, bra, ketbra, proj, bloch_vector, res, unres, max_mixed, max_entangled,
    werner_state, permutesystems
using SparseArrays

function ket(::Type{T}, val::Int, dim::Int) where T<:AbstractVector{<:Number}
    dim > 0 ? () : throw(ArgumentError("Vector dimension has to be nonnegative"))
    1 <= val <= dim ? () : throw(ArgumentError("Label have to be smaller than vector dimension"))
    ψ = T(undef, dim)
    fill!(ψ, zero(eltype(T)))
    ψ[val] = one(eltype(T))
    ψ
end

ket(::Type{<:AbstractSparseVector{T}}, val::Int, dim::Int) where T<:Number = sparsevec([val], [one(T)], dim)

ket(::Type{T}, val::Int, dim::Int) where T<:Number = ket(Vector{T}, val, dim)

"""

- `val`: non-zero entry - label.
- `dim`: length of the vector.

Return complex column vector \$|val\\rangle\$ of unit norm describing quantum state.
"""
ket(val::Int, dim::Int) = ket(ComplexF64, val, dim)


bra(::Type{T}, val::Int, dim::Int) where T<:AbstractVector{<:Number} = ket(T, val, dim)'
bra(::Type{T}, val::Int, dim::Int) where T<:Number = bra(Vector{T}, val, dim)

"""

- `val`: non-zero entry - label.
- `dim`: length of the vector

Return Hermitian conjugate \$\\langle val| = |val\\rangle^\\dagger\$ of the ket with the same label.
"""
bra(val::Int, dim::Int) = bra(ComplexF64, val, dim)

function ketbra(::Type{T}, valk::Int, valb::Int, idim::Int, odim::Int) where T<:AbstractMatrix{<:Number}
    idim > 0 && odim > 0 ? () : throw(ArgumentError("Matrix dimension has to be nonnegative"))
    1 <= valk <= idim && 1 <= valb <= odim ? () : throw(ArgumentError("Ket and bra labels have to be smaller than operator dimension"))
    ρ = T(undef, odim, idim)
    fill!(ρ, zero(eltype(T)))
    ρ[valk,valb] = one(eltype(T))
    ρ
end

ketbra(::Type{<:AbstractSparseMatrix{T}}, valk::Int, valb::Int, idim::Int, odim::Int) where T<:Number = sparse([valb], [valk], [one(T)], odim, idim)

ketbra(::Type{T}, valk::Int, valb::Int, dim::Int) where T<:AbstractMatrix{<:Number} = ketbra(T, valk, valb, dim, dim)
ketbra(::Type{T}, valk::Int, valb::Int, dim::Int) where T<:Number = ketbra(Matrix{T}, valk, valb, dim)

"""

- `valk`: non-zero entry - label.
- `valb`: non-zero entry - label.
- `dim`: length of the ket and bra vectors

# Return outer product \$|valk\\rangle\\langle vakb|\$ of states \$|valk\\rangle\$ and \$|valb\\rangle\$.
"""
ketbra(valk::Int, valb::Int, dim::Int) = ketbra(ComplexF64, valk, valb, dim)


"""
- `valk`: non-zero entry - label.
- `valb`: non-zero entry - label.
- `idim`: length of the ket vector
- `odim`: length of the bra vector

# Return outer product \$|valk\\rangle\\langle vakb|\$ of states \$|valk\\rangle\$ and \$|valb\\rangle\$.
"""
ketbra(valk::Int, valb::Int, idim::Int, odim::Int) = ketbra(Matrix{ComplexF64}, valk, valb, idim, odim)

"""

- `ket`: input column vector.

Return outer product \$|ket\\rangle\\langle ket|\$ of `ket`.
"""
proj(ψ::AbstractVector{<:Number}) = ψ * ψ'

"""

- `ρ`: input qubit density matrix.

Return the Bloch vector corresponding to the inpu quit state.
"""
function bloch_vector(ρ::AbstractMatrix{T}) where {T <: Number}
    @assert size(ρ) == (2, 2)
    x = 2real(ρ[1, 2])
    y = 2imag(ρ[1, 2])
    z = 2real(ρ[1, 1]) - 1
    T[x, y, z]
end

"""

- `ρ`: input matrix.

Returns `vec(ρ.T)`. Reshaping maps
    matrix `ρ` into a vector row by row.
"""
res(ρ::AbstractMatrix{<:Number}) = @cast x[(j, i)] := ρ[i, j] i in 1:size(ρ, 1), j in 1:size(ρ, 2)

unres(ϕ::AbstractVector{<:Number}, cols::Int) = @cast x[i, j] := ϕ[(j, i)] j in 1:cols

"""

- `ϕ`: input matrix.

Return de-reshaping of the vector into a matrix.
"""
function unres(ρ::AbstractVector{<:Number})
    dim = size(ρ, 1)
    s = isqrt(dim)
    unres(ρ, s)
end


"""

- `d`: length of the vector.

Return maximally mixed state \$\\frac{1}{d}\\sum_{i=0}^{d-1}|i\\rangle\\langle i |\$ of length \$d\$.
"""
max_mixed(d::Int) = I(d)/d

"""

- `d`: length of the vector.

Return maximally entangled state \$\\frac{1}{\\sqrt{d}}\\sum_{i=0}^{\\sqrt{d}-1}|ii\\rangle\$ of length \$\\sqrt{d}\$.
"""
function max_entangled(::Type{T}, d::Int) where T<:Number
    sd = isqrt(d)
    ρ = res(Diagonal{T}(I, sd))
    renormalize!(ρ)
    poster = convert(Vector{T}, ρ)
    poster
end

function max_entangled(::Type{<:AbstractSparseVector{T}}, d::Int) where T<:Number
    sd = isqrt(d)
    # sum |ii> for i=0..sd-1
    # |ii> -> index i*sd + i + 1 (1-based)
    # e.g. d=4, sd=2. |00>->1, |11>->4.
    indices = [i*sd + i + 1 for i in 0:sd-1]
    vals = fill(one(T), sd)
    
    vec = sparsevec(indices, vals, d)
    renormalize!(vec)
    vec
end

max_entangled(d::Int) = max_entangled(ComplexF64, d)

"""

- `d`: length of the vector.
- `α`: real number from [0, 1].

Returns [Werner state](http://en.wikipedia.org/wiki/Werner_state) given by
\$\\frac{\\alpha}{d}\\left(\\sum_{i=0}^{\\sqrt{d}-1}|ii\\rangle\\right)
\\left(\\sum_{i=0}^{\\sqrt{d}-1}\\langle ii|\\right)+
\\frac{1-\\alpha}{d}\\sum_{i=0}^{d-1}|i\\rangle\\langle i|\$.
"""
function werner_state(d::Int, α::Number)
    α > 1 || α < 0 ? throw(ArgumentError("α must be in [0, 1]")) : ()
    α * proj(max_entangled(complex(typeof(α)), d)) + (1 - α) * max_mixed(d)
end

"""

- `ρ`: input state.
- `dims`: dimensions of registers of `ρ`.
- `systems`: permuted registers.

Returns state ρ with permuted registers denoted by `systems`.
"""
function permutesystems(ρ::AbstractMatrix{T}, dims::Vector{Int}, systems::Vector{Int}) where T<:Number
    if size(ρ,1)!=size(ρ,2)
        throw(ArgumentError("Non square matrix passed to ptrace"))
    end
    if prod(dims)!=size(ρ,1)
        throw(ArgumentError("Product of dimensions does not match the shape of matrix."))
    end
    if maximum(systems) > length(dims) || minimum(systems) < 1
        throw(ArgumentError("System index out of range"))
    end
    offset = length(dims)
    perm_1 = systems
    perm_2 = [p + offset for p in perm_1]
    perm = [perm_1 ; perm_2] # vcat(perm_1 ; perm_2)
    reversed_indices = tuple(collect(length(perm):-1:1)...)
    reversed_dims = reverse(dims)
    tensor = reshape(ρ, tuple([reversed_dims ; reversed_dims]...))

    # reversed_tensor is introduced because of differences how arrays are stored and reshaped in julia and numpy
    reversed_tensor = permutedims(tensor, reversed_indices)
    reversed_transposed_tensor = permutedims(reversed_tensor, perm)
    transposed_tensor = permutedims(reversed_transposed_tensor, reversed_indices)
    return reshape(transposed_tensor, size(ρ))
end

function permutesystems(ρ::AbstractSparseMatrix{T}, dims::Vector{Int}, systems::Vector{Int}) where T<:Number
    if size(ρ,1) != size(ρ,2)
        throw(ArgumentError("Non square matrix passed to permutesystems"))
    end
    if prod(dims) != size(ρ,1)
        throw(ArgumentError("Product of dimensions does not match the shape of matrix."))
    end
    N = length(dims)
    if maximum(systems) > N || minimum(systems) < 1
        throw(ArgumentError("System index out of range"))
    end

    rev_dims = reverse(dims)
    
    I, J, V = findnz(ρ)
    new_I = Vector{Int}(undef, length(I))
    new_J = Vector{Int}(undef, length(J))
    
    # Temporary buffers
    sys_indices = Vector{Int}(undef, N)
    
    for k in 1:length(V)
        for (idx_arr, val_ptr) in ((new_I, I[k]), (new_J, J[k]))
            val = val_ptr - 1
            for step in 1:N
                d = rev_dims[step]
                rem = val % d
                val = div(val, d)
                sys_id = N - step + 1
                sys_indices[sys_id] = rem
            end
            
            new_val = 0
            current_stride = 1
            
            new_dims = dims[systems]
            rev_new_dims = reverse(new_dims)
            
            for step in 1:N
                target_sys = systems[N - step + 1]
                idx_val = sys_indices[target_sys]
                
                new_val += idx_val * current_stride
                current_stride *= rev_new_dims[step]
            end
            
            idx_arr[k] = new_val + 1
        end
    end
    
    sparse(new_I, new_J, V, size(ρ, 1), size(ρ, 2))
end

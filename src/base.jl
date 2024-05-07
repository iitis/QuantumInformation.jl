export ket, bra, ketbra, proj, bloch_vector, res, unres, max_mixed, max_entangled,
    werner_state, permutesystems

function ket(::Type{T}, val::Int, dim::Int) where T<:AbstractVector{<:Number}
    dim > 0 ? () : throw(ArgumentError("Vector dimension has to be nonnegative"))
    1 <= val <= dim ? () : throw(ArgumentError("Label have to be smaller than vector dimension"))
    ψ = T(undef, dim)
    fill!(ψ, zero(eltype(T)))
    ψ[val] = one(eltype(T))
    ψ
end

ket(::Type{T}, val::Int, dim::Int) where T<:Number = ket(Vector{T}, val, dim)

"""
$(SIGNATURES)
- `val`: non-zero entry - label.
- `dim`: length of the vector.

Return complex column vector \$|val\\rangle\$ of unit norm describing quantum state.
"""
ket(val::Int, dim::Int) = ket(ComplexF64, val, dim)


bra(::Type{T}, val::Int, dim::Int) where T<:AbstractVector{<:Number} = ket(T, val, dim)'
bra(::Type{T}, val::Int, dim::Int) where T<:Number = bra(Vector{T}, val, dim)

"""
$(SIGNATURES)
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

ketbra(::Type{T}, valk::Int, valb::Int, dim::Int) where T<:AbstractMatrix{<:Number} = ketbra(T, valk, valb, dim, dim)
ketbra(::Type{T}, valk::Int, valb::Int, dim::Int) where T<:Number = ketbra(Matrix{T}, valk, valb, dim)

"""
$(SIGNATURES)
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
$(SIGNATURES)
- `ket`: input column vector.

Return outer product \$|ket\\rangle\\langle ket|\$ of `ket`.
"""
proj(ψ::AbstractVector{<:Number}) = ψ * ψ'

"""
$(SIGNATURES)
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
$(SIGNATURES)
- `ρ`: input matrix.

Returns `vec(ρ.T)`. Reshaping maps
    matrix `ρ` into a vector row by row.
"""
res(ρ::AbstractMatrix{<:Number}) = @cast x[(j, i)] := ρ[i, j]

unres(ϕ::AbstractVector{<:Number}, cols::Int) = @cast x[i, j] := ϕ[(j, i)] j in 1:cols

"""
$(SIGNATURES)
- `ϕ`: input matrix.

Return de-reshaping of the vector into a matrix.
"""
function unres(ρ::AbstractVector{<:Number})
    dim = size(ρ, 1)
    s = isqrt(dim)
    unres(ρ, s)
end


"""
$(SIGNATURES)
- `d`: length of the vector.

Return maximally mixed state \$\\frac{1}{d}\\sum_{i=0}^{d-1}|i\\rangle\\langle i |\$ of length \$d\$.
"""
max_mixed(d::Int) = I(d)/d

"""
$(SIGNATURES)
- `d`: length of the vector.

Return maximally entangled state \$\\frac{1}{\\sqrt{d}}\\sum_{i=0}^{\\sqrt{d}-1}|ii\\rangle\$ of length \$\\sqrt{d}\$.
"""
function max_entangled(d::Int)
    sd = isqrt(d)
    ρ = res(Diagonal{ComplexF64}(I, sd))
    renormalize!(ρ)
    ρ
end

"""
$(SIGNATURES)
- `d`: length of the vector.
- `α`: real number from [0, 1].

Returns [Werner state](http://en.wikipedia.org/wiki/Werner_state) given by
\$\\frac{\\alpha}{d}\\left(\\sum_{i=0}^{\\sqrt{d}-1}|ii\\rangle\\right)
\\left(\\sum_{i=0}^{\\sqrt{d}-1}\\langle ii|\\right)+
\\frac{1-\\alpha}{d}\\sum_{i=0}^{d-1}|i\\rangle\\langle i|\$.
"""
function werner_state(d::Int, α::Float64)
    α > 1 || α < 0 ? throw(ArgumentError("α must be in [0, 1]")) : ()
    α * proj(max_entangled(d)) + (1 - α) * max_mixed(d)
end

"""
$(SIGNATURES)
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

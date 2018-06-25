function ket(::Type{Tv}, val::Int, dim::Int) where Tv<:AbstractVector{T} where T<:Number
    dim > 0 ? () : throw(ArgumentError("Vector dimension has to be nonnegative"))
    val < dim ? () : throw(ArgumentError("Label have to be smaller than vector dimmension"))
    ψ = zeros(T, dim)
    ψ[val+1] = one(T)
    ψ
end

"""
$(SIGNATURES)
- `val`: non-zero entry - label.
- `dim`: length of the vector.

Return complex column vector \$|val\\rangle\$ of unit norm describing quantum state.
"""
ket(val::Int, dim::Int) = ket(Vector{ComplexF64}, val, dim)


bra(::Type{Tv}, val::Int, dim::Int) where Tv<:AbstractVector{T} where T<:Number = ket(Tv, val, dim)'

"""
$(SIGNATURES)
- `val`: non-zero entry - label.
- `dim`: length of the vector

Return Hermitian conjugate \$\\langle val| = |val\\rangle^\\dagger\$ of the ket with the same label.
"""
bra(val::Int, dim::Int) = bra(Vector{ComplexF64}, val, dim)

function ketbra(::Type{Tm}, valk::Int, valb::Int, dim::Int) where Tm<:AbstractMatrix{T} where T<:Number
    dim > 0 ? () : throw(ArgumentError("Vector dimension has to be nonnegative"))
    valk < dim && valb < dim ? () : throw(ArgumentError("Ket and bra labels have to be smaller than operator dimmension"))
    ρ = zeros(T, dim, dim)
    ρ[valk+1,valb+1] = one(T)
    ρ
end

"""
$(SIGNATURES)
- `valk`: non-zero entry - label.
- `valb`: non-zero entry - label.
- `dim`: length of the vector

Return outer product \$\|valk\\rangle\\langle vakb|\$ of states \$\|valk\\rangle\$ and \$\|valb\\rangle\$.
"""
ketbra(valk::Int, valb::Int, dim::Int) = ketbra(Matrix{ComplexF64}, valk, valb, dim)

"""
$(SIGNATURES)
- `ket`: input column vector.

Return outer product \$|ket\\rangle\\langle ket|\$ of `ket`.
"""
proj(ψ::AbstractVector{<:Number}) = ψ * ψ'

"""
$(SIGNATURES)
- `dim`: length of the matrix.

Returns elementary matrices of dimension `dim` x `dim`.
"""
# TODO: allow rectangular matrices
base_matrices(::Type{Tm}, dim::Int) where Tm<:AbstractMatrix{<:Number} = Channel() do c
    dim > 0 ? () : error("Operator dimension has to be nonnegative")
    for i=0:dim-1, j=0:dim-1
        push!(c, ketbra(Tm, j, i, dim))
    end
end

base_matrices(dim::Int) = base_matrices(Matrix{ComplexF64}, dim)

"""
$(SIGNATURES)
- `ρ`: input matrix.

Returns `vec(ρ.T)`. Reshaping maps
    matrix `ρ` into a vector row by row.
"""
res(ρ::AbstractMatrix{<:Number}) = vec(transpose(ρ))

function unres(ϕ::AbstractVector{<:Number}, cols::Int)
    dim = length(ϕ)
    rows = div(dim, cols)
    rows*cols == dim ? () : throw(ArgumentError("Wrong number of columns"))
    transpose(reshape(ϕ, cols, rows))
end

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
max_mixed(d::Int) = eye(ComplexF64, d, d)/d

"""
$(SIGNATURES)
- `d`: length of the vector.

Return maximally entangled state \$\\frac{1}{\\sqrt{d}}\\sum_{i=0}^{\\sqrt{d}-1}|ii\\rangle\$ of length \$\\sqrt{d}\$.
"""
function max_entangled(d::Int)
    sd = isqrt(d)
    ρ = res(eye(ComplexF64, sd, sd))
    renormalize!(ρ)
    ρ
end

"""
- `d`: length of the vector.
- `α`: real number from [0, 1].

Returns [Werner state](http://en.wikipedia.org/wiki/Werner_state) given by
\$ \\frac{\\alpha}{d}\\Big(\\sum_{i=0}^{\\sqrt{d}-1}|ii\\rangle\\Big) \\Big(\\sum_{i=0}^{\\sqrt{d}-1}\\langle ii|\\Big)+ \\frac{1-\\alpha}{d}\\sum_{i=0}^{d-1}|i\\rangle\\langle i |\$.
"""
function werner_state(d::Int, α::Float64)
    α > 1 || α < 0 ? throw(ArgumentError("α must be in [0, 1]")) : ()
    α * proj(max_entangled(d)) + (1 - α) * max_mixed(d)
end


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
    tensor = reshape(ρ, tuple([dims ; dims]...))
    transposed_tensor = permutedims(tensor, perm)
    return reshape(transposed_tensor, size(ρ))
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

=#

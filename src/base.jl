function ket(::Type{Tv}, val::Int, dim::Int) where Tv<:AbstractVector{T} where T<:Number
    dim > 0 ? () : throw(ArgumentError("Vector dimension has to be nonnegative"))
    val < dim ? () : throw(ArgumentError("Label have to be smaller than vector dimmension"))
    ϕ = zeros(T, dim)
    ϕ[val+1] = one(T)
    ϕ
end

function ket(::Type{Tv}, val::Int, dim::Int) where Tv<:AbstractSparseVector{T} where T<:Number
    dim > 0 ? () : throw(ArgumentError("Vector dimension has to be nonnegative"))
    val < dim ? () : throw(ArgumentError("Label have to be smaller than vector dimmension"))
    ϕ = spzeros(T, dim)
    ϕ[val+1] = one(T)
    ϕ
end

"""
$(SIGNATURES)
- `val`: non-zero entry - label.
- `dim`: length of the vector.
- `sparse` : sparse\/dense option. Optional `sparse=false`.

Return complex column vector \$|val\\rangle\$ of unit norm describing quantum state.
"""
function ket(val::Int, dim::Int; sparse=false)
    if sparse
        return ket(SparseVector{ComplexF64}, val, dim)
    else
        return ket(Vector{ComplexF64}, val, dim)
    end
end


bra(::Type{Tv}, val::Int, dim::Int) where Tv<:AbstractVector{T} where T<:Number = ket(Tv, val, dim)'

"""
$(SIGNATURES)
- `val`: non-zero entry - label.
- `dim`: length of the vector
- `sparse` : sparse\/dense option. Optional `sparse=false`.

Return Hermitian conjugate \$\\langle val| = |val\\rangle^\\dagger\$ of the ket with the same label.
"""
function bra(val::Int, dim::Int; sparse=false)
    if sparse
        return bra(SparseVector{ComplexF64}, val, dim)
    else
        return bra(Vector{ComplexF64}, val, dim)
    end
end

function ketbra(::Type{Tm}, valk::Int, valb::Int, dim::Int) where Tm<:AbstractMatrix{T} where T<:Number
    dim > 0 ? () : throw(ArgumentError("Vector dimension has to be nonnegative"))
    valk < dim && valb < dim ? () : throw(ArgumentError("Ket and bra labels have to be smaller than operator dimmension"))
    ϕψ = zeros(T, dim, dim)
    ϕψ[valk+1,valb+1] = one(T)
    ϕψ
end

function ketbra(::Type{Tm}, valk::Int, valb::Int, dim::Int) where Tm<:AbstractSparseMatrix{T} where T<:Number
    dim > 0 ? () : throw(ArgumentError("Vector dimension has to be nonnegative"))
    valk < dim && valb < dim ? () : throw(ArgumentError("Ket and bra labels have to be smaller than operator dimmension"))
    ϕψ = spzeros(T, dim, dim)
    ϕψ[valk+1,valb+1] = one(T)
    ϕψ
end

"""
$(SIGNATURES)
- `valk`: non-zero entry - label.
- `valb`: non-zero entry - label.
- `dim`: length of the vector
- `sparse` : sparse\/dense option. Optional `sparse=false`.

Return outer product \$\|valk\\rangle\\langle vakb|\$ of states \$\|valk\\rangle\$ and \$\|valb\\rangle\$.
"""
function ketbra(valk::Int, valb::Int, dim::Int; sparse=false)
    if sparse
        return ketbra(SparseMatrixCSC{ComplexF64}, valk, valb, dim)
    else
        return ketbra(Matrix{ComplexF64}, valk, valb, dim)
    end
end

"""
$(SIGNATURES)
- `ket`: input column vector.

Return outer product \$|ket\\rangle\\langle ket|\$ of `ket`.
"""
proj(ket::AbstractVector{<:Number}) = ket * ket'

# function base_matrices(dim)
#     function _it()
#         for i=0:dim-1, j=0:dim-1
#             produce(ketbra(j, i, dim))
#         end
#     end
#     Task(_it)
# end

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


"""
$(SIGNATURES)
- `ϕ`: input matrix.

Return de-reshaping of the vector into a matrix.
"""
function unres(ϕ::AbstractVector{<:Number}, cols::Int)
    dim = length(ϕ)
    rows = div(dim, cols)
    rows*cols == dim ? () : error("Wrong number of columns")
    transpose(reshape(ϕ, cols, rows))
end

function unres(ϕ::AbstractVector{<:Number})
    dim = size(ϕ, 1)
    s = isqrt(dim)
    unres(ϕ, s)
end

"""
$(SIGNATURES)
- `kraus_list`: list of vectors.
- `ρ`: input matrix.

Return mapping of `kraus_list` on `ρ`. Krauss representation of quantum channel
\$\\Phi\$ is a set \$\\{K_i\\}_{i\\in I}\$ of bounded operators on \$\\mathcal{H}\$
such that \$\\sum_{i\\in I} K_i^\\dagger K_i = \\mathcal{1}\$.
Then \$\\Phi(\\rho)=\\sum_{i\\in I} K_i \\rho K_i^\\dagger\$.
"""
# TODO: allow different type of Kraus operators and the quantum state
function apply_kraus(kraus_list::Vector{T}, ρ::T) where {T<:AbstractMatrix{T1}} where {T1<:Number}
    # TODO: chceck if all Kraus operators are the same shape and fit the input state
    sum(k-> k*ρ*k', kraus_list)
end

"""
$(SIGNATURES)
- `d`: length of the vector.
- `sparse` : sparse\/dense option. Optional `sparse=false`.

Return maximally mixed state \$\\frac{1}{d}\\sum_{i=0}^{d-1}|i\\rangle\\langle i |\$ of length \$d\$.
"""
max_mixed(d::Int; sparse=false) = sparse ? speye(ComplexF64, d, d)/d : eye(ComplexF64, d, d)/d

"""
$(SIGNATURES)
- `d`: length of the vector.
- `sparse` : sparse\/dense option. Optional `sparse=false`.

Return maximally entangled state \$\\frac{1}{\\sqrt{d}}\\sum_{i=0}^{\\sqrt{d}-1}|ii\\rangle\$ of length \$\\sqrt{d}\$.
"""
function max_entangled(d::Int; sparse=false)
    sd = isqrt(d)
    ϕ = sparse ? res(speye(ComplexF64, sd, sd)) : res(eye(ComplexF64, sd, sd))
    renormalize!(ϕ)
    ϕ
end

"""
- `d`: length of the vector.
- `α`: real number from [0, 1].

Returns [Werner state](http://en.wikipedia.org/wiki/Werner_state) given by
\$ \\frac{\\alpha}{d}\\Big(\\sum_{i=0}^{\\sqrt{d}-1}|ii\\rangle\\Big) \\Big(\\sum_{i=0}^{\\sqrt{d}-1}\\langle ii|\\Big)+ \\frac{1-\\alpha}{d}\\sum_{i=0}^{d-1}|i\\rangle\\langle i |\$.
"""
function werner_state(d::Int, α::Float64,)
    α > 1 || α < 0 ? throw(ArgumentError("α must be in [0, 1]")) : ()
    α * proj(max_entangled(d)) + (1 - α) * max_mixed(d)
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


def permute_systems(rho, dims, systemperm):
    rho = np.asarray(rho)
    dims = list(dims)
    systemperm = list(systemperm)
    if rho.shape[0] != rho.shape[1]:
        raise Exception("Non square matrix passed to ptrace")
    if np.prod(dims) != rho.shape[0]:
        raise Exception("Product of dimensions do not match shape of matrix.")
    if not ((max(systemperm) <= len(dims) or (min(systemperm) > len(dims)))):
        raise Exception("System index out of range")
    offset = len(dims)
    perm1 = systemperm
    perm2 = map(lambda x: x + offset, perm1)
    perm = perm1 + perm2
    tensor = np.array(rho).reshape(2 * dims)
    tensor = tensor.transpose(perm)
    return np.asmatrix(tensor.reshape(rho.shape))
=#

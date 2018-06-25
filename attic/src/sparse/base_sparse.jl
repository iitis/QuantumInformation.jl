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

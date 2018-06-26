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

function ket(val::Int, dim::Int; sparse=false)
    if sparse
        return ket(SparseVector{ComplexF64}, val, dim)
    else
        return ket(Vector{ComplexF64}, val, dim)
    end
end

bra(::Type{Tv}, val::Int, dim::Int) where Tv<:AbstractVector{T} where T<:Number = ket(Tv, val, dim)'

function bra(val::Int, dim::Int; sparse=false)
    if sparse
        return bra(SparseVector{ComplexF64}, val, dim)
    else
        return bra(Vector{ComplexF64}, val, dim)
    end
end

function ketbra(::Type{Tv}, valk::Int, valb::Int, dim::Int) where Tv<:AbstractMatrix{T} where T<:Number
    dim > 0 ? () : throw(ArgumentError("Vector dimension has to be nonnegative"))
    valk < dim && valb < dim ? () : throw(ArgumentError("Ket and bra labels have to be smaller than operator dimmension"))
    ϕψ = zeros(T, dim, dim)
    ϕψ[valk+1,valb+1] = one(T)
    ϕψ
end

function ketbra(::Type{Tv}, valk::Int, valb::Int, dim::Int) where Tv<:AbstractSparseMatrix{T} where T<:Number
    dim > 0 ? () : throw(ArgumentError("Vector dimension has to be nonnegative"))
    valk < dim && valb < dim ? () : throw(ArgumentError("Ket and bra labels have to be smaller than operator dimmension"))
    ϕψ = spzeros(T, dim, dim)
    ϕψ[valk+1,valb+1] = one(T)
    ϕψ
end

function ketbra(valk::Int, valb::Int, dim::Int; sparse=false)
    if sparse
        return ketbra(SparseMatrixCSC{ComplexF64}, valk, valb, dim)
    else
        return ketbra(Matrix{ComplexF64}, valk, valb, dim)
    end
end

proj(ket::AbstractVector{T}) where T<:Number = ket * ket'

# function base_matrices(dim)
#     function _it()
#         for i=0:dim-1, j=0:dim-1
#             produce(ketbra(j, i, dim))
#         end
#     end
#     Task(_it)
# end

base_matrices(dim) = Channel() do c
    dim > 0 ? () : error("Operator dimension has to be nonnegative")
    for i=0:dim-1, j=0:dim-1
        push!(c, ketbra(j, i, dim))
    end
end

res(ρ::AbstractMatrix{T}) where T<:Number = vec(transpose(ρ))

function unres(ϕ::AbstractVector{T}, cols::Int) where T<:Number
    dim = length(ϕ)
    rows = div(dim, cols)
    rows*cols == length(ϕ) ? () : error("Wrong number of columns")
    transpose(reshape(ϕ, cols, rows))
end

function unres(ϕ::AbstractVector{T}) where T<:Number
    dim = size(ϕ, 1)
    s = isqrt(dim)
    unres(ϕ, s)
end

unres(ϕ::AbstractVector{T}, m::Int, n::Int) where T<:Number = transpose(reshape(ϕ, n, m))

# TODO: allow different type of Kraus operators and the quantum state
function apply_kraus(kraus_list::Vector{T}, ρ::T) where {T<:AbstractMatrix{T1}} where {T1<:Number}
    # TODO: chceck if all Kraus operators are the same shape and fit the input state
    sum(k-> k*ρ*k', kraus_list)
end

max_mixed(d::Int; sparse=false) = sparse ? speye(ComplexF64, d, d)/d : eye(ComplexF64, d, d)/d

function max_entangled(d::Int; sparse=false)
    sd = isqrt(d)
    ϕ = sparse ? res(speye(ComplexF64, sd, sd)) : res(eye(ComplexF64, sd, sd))
    renormalize!(ϕ)
    ϕ
end

"""
http://en.wikipedia.org/wiki/Werner_state
"""
function werner_state(d::Int, α::Float64,)
    α > 1 || α < 0 ? throw(ArgumentError("α must be in [0, 1]")) : ()
    α * proj(max_entangled(d)) + (1 - α) * max_mixed(d)
end


function permute_systems(ρ::AbstractMatrix{T}, dims::Vector{Int}, systems::Vector{Int}) where T<:Number
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
    reversed_indices = (length(perm):-1:1...)
    tensor = reshape(ρ, tuple([dims ; dims]...))
    # reversed_tensor is introduced because of differences how arrays are stored and reshaped in julia and numpy
    reversed_tensor = permutedims(tensor, reversed_indices)
    reversed_transposed_tensor = permutedims(reversed_tensor, perm)
    transposed_tensor = permutedims(reversed_transposed_tensor, reversed_indices)
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

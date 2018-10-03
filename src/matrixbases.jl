abstract type AbstractMatrixBase{T<:AbstractMatrix{<:Number}} end
struct ElementaryMatrixBase{T}<:AbstractMatrixBase{T}
    channel::Channel{T}
    dim::Int
end

function ElementaryMatrixBase{T}(::Type{T}, dim::Int)
    ElementaryMatrixBase(base_matrices(::Type{T}, dim::Int), dim)
end

# TODO: allow rectangular matrices
base_matrices(::Type{T}, dim::Int) where T<:Number = Channel() do c
    dim > 0 ? () : error("Operator dimension has to be nonnegative")
    for i=1:dim, j=1:dim
        push!(c, ketbra(T, j, i, dim))
    end
end

"""
$(SIGNATURES)
- `dim`: length of the matrix.

Returns elementary matrices of dimension `dim` x `dim`.
"""
base_matrices(dim::Int) = base_matrices(ComplexF64, dim)

"""
$(SIGNATURES)
- `dim`: dimensions of registers of `ρ`.

Returns elementary hermitian matrices of dimension dim x dim.
"""
base_hermitian_matrices(dim) = Channel(ctype=Matrix{ComplexF64}) do bhm
    for (a, b) in Base.product(0:dim-1, 0:dim-1)
        if a > b
            x = 1 / sqrt(2) * (1im * ketbra(a, b, dim) - 1im * ketbra(b, a, dim))
            push!(bhm, x)
        elseif a < b
            x = 1 / sqrt(2) * (ketbra(a, b, dim) + ketbra(b, a, dim))
            push!(bhm, x)
        else
            x = ketbra(a, b, dim)
            push!(bhm, x)
        end
    end
    bhm
end

base_generlized_pauli_matrices(d) = Channel(ctype=Matrix{ComplexF64}) do gm
    E(i,j,d) = ketbra(j,i,d)

    sm = sum([E(i,i,d) for i in 0:d-1])
    push!(gm, sm)
    for j in 0:d-2
        for k in j+1:d-1
            sm = E(j,k,d) + E(k,j,d)
            push!(gm, sm)
        end
    end

    for j in 0:d-2
        for k in j+1:d-1
            sm = -1im*(E(j,k,d) - E(k,j,d))
            push!(gm, sm)
        end
    end

    for l in 1:d-1
        sm = sqrt(2.0/(l*(l+1)))*(sum([E(j-1,j-1,d) for j in 1:l]) - l*E(l,l,d))
        push!(gm, sm)
    end
    gm
end

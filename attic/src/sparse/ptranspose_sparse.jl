function ptranspose(ρ::AbstractSparseMatrix{<:Number}, dims::Vector{Int}, sys::Int)

    if size(ρ,1)!=size(ρ,2)
        throw(ArgumentError("Non square matrix passed to ptrace"))
    end
    if prod(dims)!=size(ρ,1)
        throw(ArgumentError("Product of dimensions do not match shape of matrix."))
    end

    if sys != 1 && sys != 2
        throw(ArgumentError("sys must be 1 or 2, not $sys"))
    end

    I, J, V = findnz(ρ)
    newI = zeros(I)
    newJ = zeros(J)
    for k=1:length(I)
        i, j = number2mixedradix(I[k]-1, dims), number2mixedradix(J[k]-1, dims)
        i[sys], j[sys] = j[sys], i[sys]
        newI[k], newJ[k] = mixedradix2number(i, dims), mixedradix2number(j, dims)
    end
    sparse(newI+1, newJ+1, V, size(ρ)...)
end

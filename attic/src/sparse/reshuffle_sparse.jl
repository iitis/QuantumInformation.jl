function reshuffle(ρ::AbstractSparseMatrix{<:Number}, dims::Matrix{Int})
    dimsI =dims[1,:]
    dimsJ =dims[2,:]
    newdimsI =[dims[1, 1], dims[2, 1]]
    newdimsJ =[dims[1, 2], dims[2, 2]]
    I, J, V = findnz(ρ)
    newI = zero(I)
    newJ = zero(J)
    for k=1:length(I)
        i, j = number2mixedradix(I[k]-1, dimsI), number2mixedradix(J[k]-1, dimsJ)
        i[1], i[2], j[1], j[2] = i[1], j[1], i[2], j[2]
        newI[k], newJ[k] = mixedradix2number(i, newdimsI), mixedradix2number(j, newdimsJ)
    end
    sparse(newI+1, newJ+1, V, prod(newdimsI), prod(newdimsJ))
end

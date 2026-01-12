export reshuffle
using SparseArrays

# """
#   Performs reshuffling of indices of a matrix.
#   Given multiindexed matrix M_{(m,μ),(n,ν)} it returns
#   matrix M_{(m,n),(μ,ν)}.
# """
function reshuffle(ρ::AbstractMatrix, dims::Matrix{Int})
    m, n, μ, ν = dims
    tensor = reshape(ρ, μ, m, ν, n)
    perm = [4, 2, 3, 1]
    tensor = permutedims(tensor, perm)
    return reshape(tensor, m*n, μ*ν)
end

function reshuffle(ρ::AbstractSparseMatrix, dims::Matrix{Int})
    m, n, μ, ν = dims

    I, J, V = findnz(ρ)

    new_I = Int[]
    new_J = Int[]
    new_V = copy(V) # Values separate? No, values stay same, positions move.

    for k in 1:length(V)
        r = I[k] - 1
        c = J[k] - 1

        i_mu = r % μ
        i_m = r ÷ μ

        i_nu = c % ν
        i_n = c ÷ ν

        new_r = i_n + n * i_m
        new_c = i_nu + ν * i_mu

        push!(new_I, new_r + 1)
        push!(new_J, new_c + 1)
    end

    return sparse(new_I, new_J, V, m*n, μ*ν)
end

"""
  - `ρ`: reshuffled matrix.
    Performs reshuffling of indices of a matrix.
    Given multiindexed matrix \$M_{(m,μ),(n,ν)}\$ it returns
    matrix \$M_{(m,n),(μ,ν)}\$.
"""
function reshuffle(ρ::AbstractMatrix)
    (r, c) = size(ρ)
    sqrtr = isqrt(r)
    sqrtc = isqrt(c)
    return reshuffle(ρ, [sqrtr sqrtr; sqrtc sqrtc])
end

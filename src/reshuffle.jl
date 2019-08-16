export reshuffle

# """
#   Performs reshuffling of indices of a matrix.
#   Given multiindexed matrix M_{(m,μ),(n,ν)} it returns
#   matrix M_{(m,n),(μ,ν)}.
# """
function reshuffle(ρ::AbstractMatrix{<:Number}, dims::Matrix{Int})
  m, n, μ, ν  = dims
  tensor = reshape(ρ, μ, m, ν, n)
  perm = [4, 2, 3, 1]
  tensor = permutedims(tensor,  perm)
  reshape(tensor, m*n, μ*ν)
end

"""
  $(SIGNATURES)
  - `ρ`: reshuffled matrix.
  Performs reshuffling of indices of a matrix.
  Given multiindexed matrix \$M_{(m,μ),(n,ν)}\$ it returns
  matrix \$M_{(m,n),(μ,ν)}\$.
"""
function reshuffle(ρ::AbstractMatrix{<:Number})
    (r, c) = size(ρ)
    sqrtr = isqrt(r)
    sqrtc = isqrt(c)
    reshuffle(ρ, [sqrtr sqrtr; sqrtc sqrtc])
end

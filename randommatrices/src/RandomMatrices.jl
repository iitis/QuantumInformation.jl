module RandomMatrices
using LinearAlgebra
import Base: rand, size
using Random: GLOBAL_RNG, AbstractRNG

export rand, size, QIContinuousMatrixDistribution

abstract type QIContinuousMatrixDistribution; end

rand(c::QIContinuousMatrixDistribution) = rand(GLOBAL_RNG, c)

include("ginibre.jl")
include("circular.jl")
include("wigner.jl")
include("wishart.jl")
end

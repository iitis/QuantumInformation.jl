module RandomMatrices
using LinearAlgebra
import Base: rand
import Distributions: ContinuousMatrixDistribution

export rand

include("ginibre.jl")
include("circular.jl")
include("wigner.jl")
include("wishart.jl")
end

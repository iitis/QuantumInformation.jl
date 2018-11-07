module RandomMatrices
using LinearAlgebra
import Base: rand
import Distributions: ContinuousMatrixDistribution

export CircularEnsemble, COE, CUE, CSE, CircularRealEnsemble,
CircularQuaternionEnsemble,
GinibreEnsemble,
WishartEnsemble,
WignerEnsemble,
rand

include("ginibre.jl")
include("circular.jl")
include("wigner.jl")
include("wishart.jl")
end

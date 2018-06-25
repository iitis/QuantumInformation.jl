module RandomMatrices
import Base: rand
import Distributions: ContinuousMatrixDistribution

export CircularEnsemble, COE, CUE, CSE, CircularRealEnsemble,
CircularQuaternionEnsemble,
GinibreEnsemble,
WishartEnsemble,
WignerEnsemble,
rand

include("circular.jl")
include("ginibre.jl")
include("wigner.jl")
include("wishart.jl")
end

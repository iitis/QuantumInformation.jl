"""
Main module for `QuantumInformation.jl` -- a Julia package for numerical computation in quantum information theory.
"""
module QuantumInformation
using LinearAlgebra
using DocStringExtensions
using TensorOperations
using Convex, SCS
using Random: AbstractRNG, GLOBAL_RNG

import Base: convert, size, length, kron, *, rand, show

const ⊗ = kron

export ⊗

include("../randommatrices/src/RandomMatrices.jl")
using .RandomMatrices
eval(Expr(:export, names(RandomMatrices)...))

include("base.jl")
include("randomqobjects.jl")
include("gates.jl")
include("utils.jl")
include("channels.jl")
include("functionals.jl")
include("reshuffle.jl")
include("ptrace.jl")
include("ptranspose.jl")
include("convex.jl")
include("matrixbases.jl")

end # module

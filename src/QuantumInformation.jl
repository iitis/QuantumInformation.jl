"""
Main module for `QuantumInformation.jl` -- a Julia package for numerical computation in quantum information theory.
"""
module QuantumInformation
using LinearAlgebra
using DocStringExtensions
using TensorOperations
using TensorCast
using Convex, SCS
using Random: AbstractRNG, GLOBAL_RNG

import Base: convert, size, length, kron, *, show, rand
import LinearAlgebra: I

const ⊗ = kron

export ⊗

using MatrixEnsembles
eval(Expr(:export, names(MatrixEnsembles)...))

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

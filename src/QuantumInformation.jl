"""
Main module for `QuantumInformation.jl` -- a Julia package for numerical computation in quantum information theory.
"""
module QuantumInformation
using LinearAlgebra
using DocStringExtensions
using TensorOperations
using Convex, SCS
using Random: AbstractRNG, GLOBAL_RNG

import Base: convert, size, length, kron, *, show, rand
import LinearAlgebra: I

const ⊗ = kron

export ⊗

include("../randommatrices/src/RandomMatrices.jl")
using .RandomMatrices
eval(Expr(:export, names(RandomMatrices)...))

using Requires
@init @require CuArrays = "3a865a2d-5b23-5a0f-bc46-62713ec82fae" include("../curandommatrices/src/CuRandomMatrices.jl")
@init @require CuArrays = "3a865a2d-5b23-5a0f-bc46-62713ec82fae" using ..CuRandomMatrices
@init @require CuArrays = "3a865a2d-5b23-5a0f-bc46-62713ec82fae" export curand

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

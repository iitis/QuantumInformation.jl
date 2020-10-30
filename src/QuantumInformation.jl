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

using Requires
function __init__()
    @require CUDA="052768ef-5323-5732-b1bb-66c8b64840ba" begin
        if CUDA.functional() && CUDA.has_cutensor()
            const CuArray = CUDA.CuArray
            const CuVector = CUDA.CuVector
            const CuMatrix = CUDA.CuMatrix
            const CuSVD = CUDA.CUSOLVER.CuSVD
            const CuQR = CUDA.CUSOLVER.CuQR
            using MatrixEnsembles
            # CUDA.allowscalar(false)
            include("cuda/randomqobjects.jl") 
        end
    end
end

end # module

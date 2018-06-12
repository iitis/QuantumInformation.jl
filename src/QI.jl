"""
Main module for `QI.jl` -- a Julia package for numerical computation in quantum information theory.
"""

module QI
if VERSION<=v"0.7"
    const ComplexF64 = Complex128
else
    using LinearAlgebra
    using SparseArrays
end

using Compat, DocStringExtensions
import Compat.Markdown

const ⊗ = kron
export ket, bra, ketbra, proj, base_matrices,
res, unres,
kraus_to_superoperator, channel_to_superoperator, apply_kraus,
ptrace, ptranspose, reshuffle,
max_mixed, max_entangled, werner_state,
number2mixedradix, mixedradix2number,
trace_norm, trace_distance, hs_norm, hs_distance,
fidelity_sqrt, fidelity, gate_fidelity,
shannon_entropy, quantum_entropy, relative_entropy, kl_divergence, js_divergence,
bures_distance, bures_angle, superfidelity,
negativity, log_negativity, ppt,
diamond_norm, diamond_distance,
random_ket, random_ket!,
random_GOE, random_GUE,
random_ginibre_matrix!, random_ginibre_matrix,
random_mixed_state!, random_mixed_state_hs, random_mixed_state,
random_dynamical_matrix!, random_dynamical_matrix,
random_jamiolkowski_state!, random_jamiolkowski_state,
random_unitary, random_orthogonal, random_isometry,
funcmh, funcmh!, renormalize!, random_ball,
sx,sy,sz, qft, hadamard, grover,
⊗

include("base.jl")
include("randommatrix.jl")
include("randomstate.jl")
include("gates.jl")
include("utils.jl")
include("channels.jl")
include("functionals.jl")
include("reshuffle.jl")
include("ptrace.jl")
include("ptranspose.jl")

end # module

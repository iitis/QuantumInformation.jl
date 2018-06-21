__precompile__()
module QI
if VERSION>v"0.7.0-DEV"
    using LinearAlgebra
    using SparseArrays
else
    const ComplexF64 = Complex128
end

const ⊗ = kron
export ket, bra, ketbra, proj, base_matrices,
res, unres,
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
sx,sy,sz, qft, hadamard, grover,⊗,
kraus_is_CPTP,
kraus_to_superoperator,
kraus_to_stinespring,
kraus_to_dynamical_matrix,
superoperator_to_kraus,
superoperator_to_stinespring,
superoperator_to_dynamical_matrix,
dynamical_matrix_to_kraus,
dynamical_matrix_to_stinespring,
dynamical_matrix_to_superoperator,
apply_channel_dynamical_matrix,
apply_channel_kraus,
apply_channel_superoperator,
apply_channel_stinespring,
channel_to_superoperator,
isidentity, ispositive

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

if VERSION>v"0.7.0-DEV"
    # Convex.jl does not support julia 0.7 yet
else
    include("convex.jl")
end

end # module

module QI
⊗ = kron
export ket, bra, ketbra, proj, base_matrices,
res, unres,
kraus_to_superoperator, channel_to_superoperator, apply_kraus,
ptrace, ptranspose, reshuffle,
number2mixedradix, mixedradix2number,
trace_distance, fidelity_sqrt, fidelity,
shannon_entropy, entropy,
random_ket, random_ket!,
random_ginibre_matrix!, random_ginibre_matrix,
random_mixed_state_hs!, random_mixed_state_hs,
random_dynamical_matrix!, random_dynamical_matrix,
random_jamiolkowski_state!, random_jamiolkowski_state,
random_unitary,
funcmh, funcmh!, renormalize!,
sx,sy,sz,
⊗

include("base.jl")
include("randommatrix.jl")
include("randomstate.jl")
include("gates.jl")
include("utils.jl")

sx = Complex128[0 1; 1 0]
sy = Complex128[0 1im; -1im 0]
sz = Complex128[1 0; 0 -1]

end # module

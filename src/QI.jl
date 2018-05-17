module QI
if VERSION<=v"0.7"
    const ComplexF64 = Complex128
else
    using LinearAlgebra
    using SparseArrays
end

const ⊗ = kron
export ket, bra, ketbra, proj, base_matrices,
res, unres,
kraus_to_superoperator, channel_to_superoperator, apply_kraus,
ptrace, ptranspose, reshuffle,
number2mixedradix, mixedradix2number,
trace_distance, fidelity_sqrt, fidelity,
shannon_entropy, entropy,
random_ket, random_ket!,
random_GOE, random_GUE,
random_ginibre_matrix!, random_ginibre_matrix,
random_mixed_state!, random_mixed_state_hs, random_mixed_state,
random_dynamical_matrix!, random_dynamical_matrix,
random_jamiolkowski_state!, random_jamiolkowski_state,
random_unitary, random_orthogonal,
funcmh, funcmh!, renormalize!,
sx,sy,sz,
⊗

include("base.jl")
include("randommatrix.jl")
include("randomstate.jl")
include("gates.jl")
include("utils.jl")
include("channels.jl")
include("functionals.jl")

end # module

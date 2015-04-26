module QI

export ket, bra, ketbra, proj, base_matrices,
res, unres,
kraus_to_superoperator, channel_to_superoperator, apply_kraus,
ptrace, reshuffle,
trace_distance, fidelity_sqrt, fidelity,
shannon_entropy, entropy,
random_ket

include("base.jl")
include("randommatrix.jl")
include("randomstate.jl")
include("gates.jl")

sx = Complex128[0 1; 1 0]
sy = Complex128[0 1im; -1im 0]
sz = Complex128[1 0; 0 -1]

end # module

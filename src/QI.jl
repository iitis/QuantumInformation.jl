module QI

include("base.jl")
include("randommatrix.jl")
include("randomstate.jl")
include("gates.jl")

sx = Complex128[0 1; 1 0]
sy = Complex128[0 1im; -1im 0]
sz = Complex128[1 0; 0 -1]

end # module

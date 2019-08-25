################################################################################
# Channels permutations
################################################################################
# function permutesystems(Φ::T, idims::Vector{Int}, odims::Vector{Int}, perm::Vector{Int}) where {T<:AbstractQuantumOperation{TM}, TM<:AbstractMatrix{<:Number}}
#     error("Not implemented")
#     ko = KrausOperators(Φ)
#     for k in ko
#         #permutesystems(k, sdims, )
#     end
# end

include("channels/constructors.jl")
include("channels/applications.jl")
include("channels/conversions.jl")
include("channels/compositions.jl")
include("channels/predicates.jl")
include("channels/io.jl")
include("channels/misc.jl")
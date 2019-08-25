################################################################################
# Channels printing
################################################################################
function Base.show(io::IO, Φ::AbstractQuantumOperation{<:Matrix{<:Number}})
    println(io, typeof(Φ))
    println(io, "    dimensions: ($(Φ.idim), $(Φ.odim))")
    if :matrix in fieldnames(typeof(Φ))
        print(io, "    ")
        print(io, Φ.matrix)
    end
    if :matrices in fieldnames(typeof(Φ))
        for (i,m) in enumerate(Φ.matrices)
            print(io, "    ")
            print(io, m)
            if i < length(Φ.matrices) println(io, "") end
        end
    end
end
# Base.show(io::IO, m::MIME"text/plain", Φ::AbstractQuantumOperation{<:Matrix{<:Number}}) = show(io, m, Φ)
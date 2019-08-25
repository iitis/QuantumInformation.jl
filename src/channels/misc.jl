export represent

################################################################################
# size() function
################################################################################
size(Φ::AbstractQuantumOperation) = (Φ.idim, Φ.odim)


################################################################################
# represent() function
################################################################################
for qop in (:KrausOperators, :POVMMeasurement)
    @eval represent(Φ::$qop) = Φ.matrices
end

for qop in (:SuperOperator, :DynamicalMatrix, :Stinespring,
            :UnitaryChannel, :PostSelectionMeasurement)
    @eval represent(Φ::$qop) = Φ.matrix
end

represent(Φ::IdentityChannel{T}) where T<:Matrix{<:Number} = Matrix{T}(I, Φ.idim, Φ.idim)
represent(Φ::IdentityChannel) = represent(IdentityChannel{Matrix{ComplexF64}}())
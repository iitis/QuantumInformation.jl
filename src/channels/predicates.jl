export ispovm, iseffect, iscp, istp, istni, iscptp, iscptni, isidentity,
    ispositive


################################################################################
# CP
################################################################################
"""

- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is completely positive.
"""
function iscp end

function iscp(Φ::KrausOperators{<:AbstractMatrix{<:Number}}, atol=1e-13)
    # by definition Kraus operators represent a CP map
    true
end

function iscp(Φ::SuperOperator{T}, atol=1e-13) where T<:AbstractMatrix{<:Number}
    iscp(convert(DynamicalMatrix{T}, Φ), atol)
end

function iscp(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}, atol=1e-13)
    ispositive(Φ.matrix, atol)
end

function iscp(Φ::Stinespring{<:AbstractMatrix{<:Number}}, atol=1e-13)
    # by definition Stinespring operator represents a CP map(?)
    true
end

function iscp(Φ::UnitaryChannel, atol=1e-13)
    # by definition Unitary operator represents a CP map
    true
end


################################################################################
# TNI
################################################################################
"""

- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is trace non-increasing.
"""
function istni end

function istni(Φ::KrausOperators{<:AbstractMatrix{<:Number}}, atol=1e-13)
    cr = sum(k'*k for k in Φ.matrices)
    ispositive(one(cr) - cr, atol)
end

function istni(Φ::SuperOperator{T}, atol=1e-13) where T<:AbstractMatrix{<:Number}
    istni(convert(DynamicalMatrix{T}, Φ), atol)
end

function istni(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}, atol=1e-13)
    pt = ptrace(Φ.matrix, [Φ.odim, Φ.idim], [1])
    ispositive(one(pt) - pt, atol)
end

function istni(Φ::Stinespring{<:AbstractMatrix{<:Number}}, atol=1e-13)
    u = Φ.matrix
    m = u'*u
    ispositive(one(m) - m, atol)
end

function istni(Φ::UnitaryChannel, atol=1e-13)
    iscptp(Φ, atol)
end

################################################################################
# TP
################################################################################
"""

- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is trace preserving.
"""
function istp end

function istp(Φ::KrausOperators{<:AbstractMatrix{<:Number}}, atol=1e-13)
    cr = sum(k'*k for k in Φ.matrices)
    isidentity(cr, atol)
end

function istp(Φ::SuperOperator{T}, atol=1e-13) where T<:AbstractMatrix{<:Number}
    istp(convert(DynamicalMatrix{T}, Φ), atol)
end

function istp(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}, atol=1e-13)
    pt = ptrace(Φ.matrix, [Φ.odim, Φ.idim], [1])
    isidentity(pt, atol)
end

function istp(Φ::Stinespring{<:AbstractMatrix{<:Number}}, atol=1e-13)
    u = Φ.matrix
    isidentity(u'*u, atol)
end

function istp(Φ::UnitaryChannel, atol=1e-13)
    u = Φ.matrix
    isidentity(u'*u, atol) && isidentity(u*u', atol)
end

################################################################################
# CPTP, CPTNI
################################################################################
"""

- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is completely positive and trace preserving.
"""
function iscptp(Φ::AbstractQuantumOperation, atol=1e-13)
    iscp(Φ, atol) && istp(Φ, atol)
end

"""

- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is completely positive and trace non-increasing.
"""
function iscptni(Φ::AbstractQuantumOperation, atol=1e-13)
    iscp(Φ, atol) && istni(Φ, atol)
end

################################################################################
# Measurements
################################################################################
"""

- `Φ`: POVM Measurement.

Checks if a set of matrices forms a valid POVM (Positive Operator-Valued Measure).
"""
function ispovm(Φ::POVMMeasurement{<:AbstractMatrix{<:Number}})
    isidentity(sum(Φ.matrices))
end

"""

- `Φ`: Post-selection Measurement.

Checks if a matrix represents a valid quantum effect (0 <= E <= I).
"""
function iseffect(Φ::PostSelectionMeasurement{<:AbstractMatrix{<:Number}})
    e = Φ.matrix
    m = e'*e
    ispositive(one(m) - m)
end

function Base.isapprox(Φ1::AbstractQuantumOperation, Φ2::AbstractQuantumOperation; kwargs...)
    T = complex(promote_type(eltype(Φ1), eltype(Φ2)))
    s1 = convert(SuperOperator{Matrix{T}}, Φ1)
    s2 = convert(SuperOperator{Matrix{T}}, Φ2)
    return isapprox(s1.matrix, s2.matrix; kwargs...)
end
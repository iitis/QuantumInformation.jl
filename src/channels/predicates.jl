export ispovm, iseffect, iscp, istp, istni, iscptp, iscptni, isidentity,
    ispositive


################################################################################
# CP
################################################################################
"""
$(SIGNATURES)
- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is completely positive.
"""
function iscp end

function iscp(Φ::KrausOperators{<:AbstractMatrix{<:Number}}; atol=1e-13)
    # by definition Kraus operators represent a CP map
    true
end

function iscp(Φ::SuperOperator{T}; atol=1e-13) where T<:AbstractMatrix{<:Number}
    iscp(convert(DynamicalMatrix{T}, Φ), atol=atol)
end

function iscp(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}; atol=1e-13)
    ispositive(Φ.matrix)
end

function iscp(Φ::Stinespring{<:AbstractMatrix{<:Number}}; atol=1e-13)
    # by definition Stinespring operator represents a CP map(?)
    true
end

function iscp(Φ::UnitaryChannel; atol=1e-13)
    # by definition Unitary operator represents a CP map
    true
end


################################################################################
# TNI
################################################################################
"""
$(SIGNATURES)
- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is trace non-increasing.
"""
function istni end

function istni(Φ::KrausOperators{<:AbstractMatrix{<:Number}}; atol=1e-13)
    cr = sum(k'*k for k in Φ.matrices)
    ispositive(one(cr) - cr, atol=atol)
end

function istni(Φ::SuperOperator{T}; atol=1e-13) where T<:AbstractMatrix{<:Number}
    istni(convert(DynamicalMatrix{T}, Φ))
end

function istni(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}; atol=1e-13)
    pt = ptrace(Φ.matrix, [Φ.odim, Φ.idim], [1])
    ispositive(one(pt) - pt, atol=atol)
end

function istni(Φ::Stinespring{<:AbstractMatrix{<:Number}}; atol=1e-13)
    u = Φ.matrix
    m = u'*u
    ispositive(one(m) - m, atol=atol)
end

function istni(Φ::UnitaryChannel; atol=1e-13)
    iscptp(Φ, atol=atol)
end

################################################################################
# TP
################################################################################
"""
$(SIGNATURES)
- `Φ`: A subtype of AbstractQuantumOperation.
- `atol`: tolerance of approximation.

Checks if an object is trace preserving.
"""
function istp end

function istp(Φ::KrausOperators{<:AbstractMatrix{<:Number}}; atol=1e-13)
    cr = sum(k'*k for k in Φ.matrices)
    isidentity(cr, atol=atol)
end

function istp(Φ::SuperOperator{T}; atol=1e-13) where T<:AbstractMatrix{<:Number}
    istp(convert(DynamicalMatrix{T}, Φ), atol=atol)
end

function istp(Φ::DynamicalMatrix{<:AbstractMatrix{<:Number}}; atol=1e-13)
    pt = ptrace(Φ.matrix, [Φ.odim, Φ.idim], [1])
    isidentity(pt, atol=atol)
end

function istp(Φ::Stinespring{<:AbstractMatrix{<:Number}}; atol=1e-13)
    u = Φ.matrix
    isidentity(u'*u, atol=atol)
end

function istp(Φ::UnitaryChannel; atol=1e-13)
    u = Φ.matrix
    isidentity(u'*u, atol=atol) && isidentity(u*u', atol=atol)
end

################################################################################
# CPTP, CPTNI
################################################################################
function iscptp(Φ::AbstractQuantumOperation; atol=1e-13)
    iscp(Φ, atol=atol) && istp(Φ, atol=atol)
end

function iscptni(Φ::AbstractQuantumOperation; atol=1e-13)
    iscp(Φ, atol=atol) && istni(Φ, atol=atol)
end

################################################################################
# Measurements
################################################################################
function ispovm(Φ::POVMMeasurement{<:AbstractMatrix{<:Number}})
    isidentity(sum(Φ.matrices))
end

function iseffect(Φ::PostSelectionMeasurement{<:AbstractMatrix{<:Number}})
    e = Φ.matrix
    m = e'*e
    ispositive(one(m) - m)
end
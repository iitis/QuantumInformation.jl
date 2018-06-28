struct CircularEnsemble{β} <: ContinuousMatrixDistribution
    d::Int

    function CircularEnsemble{β}(d::Int) where β
        β == 4 && mod(d, 2) == 1 ? throw(ArgumentError("Dim must even")) : ()
        new(d)
    end
end

const COE = CircularEnsemble{1}
const CUE = CircularEnsemble{2}
const CSE = CircularEnsemble{4}

function _qr_fix(z::AbstractMatrix)
    q, r = qr(z)
    d = diag(r)
    ph = d./abs.(d)
    q.*repmat(ph, 1, size(q, 1))'
end

function rand(c::COE)
    z = rand(GinibreEnsemble(c.d))
    u = _qr_fix(z)
    transpose(u)*u
end

function rand(c::CUE)
    z = rand(GinibreEnsemble(c.d))
    u = _qr_fix(z)
    u
end

function rand(c::CSE)
    z = rand(GinibreEnsemble(c.d))
    u = _qr_fix(z)
    #TODO this does not require matrix multiplication
    a = diagm(-ones(c.d-1), 1) + diagm(ones(c.d-1), -1)
    ur = -a*transpose(u)*a
    ur*u
end

struct CircularRealEnsemble <: ContinuousMatrixDistribution
    d::Int
end

function rand(c::CircularRealEnsemble)
    z = rand(GinibreEnsemble{1}(c.d))
    _qr_fix(z)
end

struct CircularQuaternionEnsemble <: ContinuousMatrixDistribution
    d::Int
end

function rand(c::CircularQuaternionEnsemble)
    z = rand(GinibreEnsemble{4}(c.d))
    _qr_fix(z)
end

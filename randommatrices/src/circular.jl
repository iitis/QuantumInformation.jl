export HaarIsometry
struct CircularEnsemble{β} <: ContinuousMatrixDistribution
    d::Int
    g::GinibreEnsemble{2}

    function CircularEnsemble{β}(d::Int) where β
        β == 4 && mod(d, 2) == 1 ? throw(ArgumentError("Dim must even")) : ()
        g = GinibreEnsemble{2}(d)
        new(d, g)
    end
end

const COE = CircularEnsemble{1}
const CUE = CircularEnsemble{2}
const CSE = CircularEnsemble{4}

function _qr_fix(z::AbstractMatrix)
    q, r = qr(z)
    d = diag(r)
    ph = d./abs.(d)
    q = Matrix(q)
    idim = size(r, 1)
    for i=1:idim
        q[:, i] .*= ph[i]
    end
    q[:, 1:idim]
end

function rand(c::COE)
    z = rand(c.g)
    u = _qr_fix(z)
    ut = diagm(ones(c.d-1), 1) + diagm(-ones(c.d-1), -1)
    ut*u
end

function rand(c::CUE)
    z = rand(c.g)
    u = _qr_fix(z)
    u
end

function rand(c::CSE)
    z = rand(c.g)
    u = _qr_fix(z)
    #TODO this does not require matrix multiplication
    ur = diagm(-ones(c.d-1), 1) + diagm(ones(c.d-1), -1)
    ur*u
end

struct CircularRealEnsemble <: ContinuousMatrixDistribution
    d::Int
    g::GinibreEnsemble{1}

    function CircularRealEnsemble(d::Int)
        g = GinibreEnsemble{1}(d)
        new(d, g)
    end
end

function rand(c::CircularRealEnsemble)
    z = rand(c.g)
    _qr_fix(z)
end

struct CircularQuaternionEnsemble <: ContinuousMatrixDistribution
    d::Int
    g::GinibreEnsemble{4}

    function CircularRealEnsemble(d::Int)
        g = GinibreEnsemble{4}(d)
        new(d, g)
    end
end

function rand(c::CircularQuaternionEnsemble)
    z = rand(c.g)
    _qr_fix(z)
end


struct HaarIsometry <: ContinuousMatrixDistribution
    idim::Int
    odim::Int
    g::GinibreEnsemble{2}

    function HaarIsometry(idim::Int, odim::Int)
        idim <= odim || throw(ArgumentError("idim can't be greater than odim"))
        g = GinibreEnsemble{2}(odim, idim)
        new(idim, odim, g)
    end
end

function rand(c::HaarIsometry)
    z = rand(c.g)
    _qr_fix(z)
end
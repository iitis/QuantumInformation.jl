export CircularEnsemble, COE, CUE, CSE, CircularRealEnsemble,
    CircularQuaternionEnsemble, HaarIsometry
struct CircularEnsemble{β} <: QIContinuousMatrixDistribution
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

function _qr_fix!(z::AbstractMatrix)
    q, r = qr!(z)
    d = diag(r)
    ph = d./abs.(d)
    q = Matrix(q)
    idim = size(r, 1)
    for i=1:idim
        q[:, i] .*= ph[i]
    end
    q[:, 1:idim]
end

function _qr_fix(z::AbstractMatrix)
    a = copy(z)
    _qr_fix!(a)
end

function rand(rng::AbstractRNG, c::COE)
    z = rand(rng, c.g)
    u = _qr_fix!(z)
    transpose(u)*u
end

function rand(rng::AbstractRNG, c::CUE)
    z = rand(rng, c.g)
    u = _qr_fix!(z)
    u
end

function rand(rng::AbstractRNG, c::CSE)
    z = rand(rng, c.g)
    u = _qr_fix!(z)
    ur = cat([[0 -1; 1 0] for _=1:c.d÷2]..., dims=[1,2])
    ur*u*ur'*transpose(u)
end

struct CircularRealEnsemble <: QIContinuousMatrixDistribution
    d::Int
    g::GinibreEnsemble{1}

    function CircularRealEnsemble(d::Int)
        g = GinibreEnsemble{1}(d)
        new(d, g)
    end
end

function rand(rng::AbstractRNG, c::CircularRealEnsemble)
    z = rand(rng, c.g)
    _qr_fix!(z)
end

struct CircularQuaternionEnsemble <: QIContinuousMatrixDistribution
    d::Int
    g::GinibreEnsemble{4}

    function CircularQuaternionEnsemble(d::Int)
        g = GinibreEnsemble{4}(d)
        new(d, g)
    end
end

function rand(rng::AbstractRNG, c::CircularQuaternionEnsemble)
    z = rand(rng, c.g)
    _qr_fix!(z)
end


struct HaarIsometry <: QIContinuousMatrixDistribution
    idim::Int
    odim::Int
    g::GinibreEnsemble{2}

    function HaarIsometry(idim::Int, odim::Int)
        idim <= odim || throw(ArgumentError("idim can't be greater than odim"))
        g = GinibreEnsemble{2}(odim, idim)
        new(idim, odim, g)
    end
end

function rand(rng::AbstractRNG, c::HaarIsometry)
    z = rand(rng, c.g)
    _qr_fix!(z)
end
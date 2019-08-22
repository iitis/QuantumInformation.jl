export WignerEnsemble

struct WignerEnsemble{β} <: QIContinuousMatrixDistribution
    d::Int
    g::GinibreEnsemble{β}

    function WignerEnsemble{β}(d::Int) where β
        β == 4 && mod(d, 2) == 1 ? throw(ArgumentError("Dim must even")) : ()
        g = GinibreEnsemble{β}(d)
        new(d, g)
    end
end

WignerEnsemble(d::Int) = WignerEnsemble{2}(d)

function rand(rng::AbstractRNG, w::WignerEnsemble{β}) where β
    z = rand(rng, w.g)
    (z + z') / 2sqrt(2β * w.d)
end
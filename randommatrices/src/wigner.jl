struct WignerEnsemble{β} <: ContinuousMatrixDistribution
    d::Int

    function WignerEnsemble{β}(d::Int) where β
        β == 4 && mod(d, 2) == 1 ? throw(ArgumentError("Dim must even")) : ()
        new(d)
    end
end

WignerEnsemble(d::Int) = WignerEnsemble{2}(d)

function rand(w::WignerEnsemble{β}) where β
    z = rand(GinibreEnsemble{β}(w.d))
    (z + z') / 2sqrt(2β * w.d)
end

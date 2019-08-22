function curand(w::WignerEnsemble{β}) where β
    z = curand(rng, GinibreEnsemble{β}(w.d))
    (z + z') / 2sqrt(2β * w.d)
end
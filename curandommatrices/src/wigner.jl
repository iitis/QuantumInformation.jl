function curand(w::WignerEnsemble{β}) where β
    z = curand(w.g)
    (z + z') / 2sqrt(2β * w.d)
end
export WishartEnsemble

struct WishartEnsemble{β, K} <: QIContinuousMatrixDistribution
    d::Int
    g::GinibreEnsemble{β}

    function WishartEnsemble{β, K}(d::Int) where {β, K}
        n = round(Int, K*d)
        K*d == n ? () : throw(ArgumentError("K*d is not and integer"))
        g = GinibreEnsemble{β}(d, n)
        new(d, g)
    end
end

 WishartEnsemble{β}(d::Int) where β = WishartEnsemble{β, 1}(d)
 WishartEnsemble(d::Int) = WishartEnsemble{2}(d)

function rand(rng::AbstractRNG, w::WishartEnsemble{β, K}) where {β, K}
    
    z = rand(rng, w.g)/sqrt(2β * w.d)
    z*z'
end

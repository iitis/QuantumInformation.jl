struct WishartEnsemble{β, K} <: ContinuousMatrixDistribution
    d::Int

    function WishartEnsemble{β, K}(d::Int) where {β, K}
        K*d == round(Int, K*d) ? () : throw(ArgumentError("K*d is not and integer"))
        new(d)
    end
end

 WishartEnsemble{β}(d::Int) where β = WishartEnsemble{β, 1}(d)
 WishartEnsemble(d::Int) = WishartEnsemble{2}(d)

function rand(w::WishartEnsemble{β, K}) where {β, K}
    n = Int(K*w.d)
    z = rand(GinibreEnsemble{β}(w.d, n))/sqrt(2β * w.d)
    z*z'
end

export GinibreEnsemble

struct GinibreEnsemble{β} <: QIContinuousMatrixDistribution
    m::Int
    n::Int

    function GinibreEnsemble{β}(m::Int, n::Int) where β
        β == 4 && (m%2 == 1 || n%2 == 1) ? throw(ArgumentError("Dim must even")) : ()
        new(m, n)
    end
end

GinibreEnsemble{β}(d::Int) where β = GinibreEnsemble{β}(d, d)
GinibreEnsemble(m::Int, n::Int) = GinibreEnsemble{2}(m, n)
GinibreEnsemble(d::Int) = GinibreEnsemble(d, d)

rand(rng::AbstractRNG, g::GinibreEnsemble{1}) = randn(rng, g.m, g.n)
rand(rng::AbstractRNG, g::GinibreEnsemble{2}) = randn(rng, g.m, g.n)+1im*randn(rng, g.m, g.n)

function rand(rng::AbstractRNG, g::GinibreEnsemble{4})
    # TODO: fix dimensions of blocks
    q0=randn(rng, g.m, g.n)
    q1=randn(rng, g.m, g.n)
    q2=randn(rng, g.m, g.n)
    q3=randn(rng, g.m, g.n)
    [q0+1im*q1 q2+1im*q3; -q2+1im*q3 q0-1im*q1]
end

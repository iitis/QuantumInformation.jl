struct GinibreEnsemble{β} <: ContinuousMatrixDistribution
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

rand(g::GinibreEnsemble{1}) = randn(g.m, g.n)
rand(g::GinibreEnsemble{2}) = randn(g.m, g.n)+1im*randn(g.m, g.n)

function rand(g::GinibreEnsemble{4})
    # TODO: fix dimensions of blocks
    error("Not implemented")
    q0=randn(g.m, g.d)
    q1=randn(g.d, g.d)
    q2=randn(g.d, g.d)
    q3=randn(g.d, g.d)
    [q0 + 1im * q1 q2 + 1im * q3; -q2 + 1im * q3 q0 - 1im * q1]
end

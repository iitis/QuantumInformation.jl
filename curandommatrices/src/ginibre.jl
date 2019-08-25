curand(g::GinibreEnsemble{1}) = CuArrays.randn(g.m, g.n)
curand(g::GinibreEnsemble{2}) = CuArrays.randn(g.m, g.n)+1im*CuArrays.randn(g.m, g.n)

function curand(g::GinibreEnsemble{4})
    q0=CuArrays.randn(g.m, g.n)
    q1=CuArrays.randn(g.m, g.n)
    q2=CuArrays.randn(g.m, g.n)
    q3=CuArrays.randn(g.m, g.n)
    [q0+1im*q1 q2+1im*q3; -q2+1im*q3 q0-1im*q1]
end

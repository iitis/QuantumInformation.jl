import Random: seed!
@testset "Convex" begin

@testset "diamond norm" begin
    seed!(42)
    c = ChoiJamiolkowskiMatrices(3, 4)
    Φ = rand(c)
    @test norm_diamond(Φ) ≈ 1 atol=1e-6
end

@testset "diamond distance" begin
    d = 4
    U1 = qft(d)
    U2 = Matrix{ComplexF64}(I, d, d)
    J1 = cat([proj(U1[:, i]) for i=1:d]..., dims=[1, 2])
    J2 = cat([proj(U2[:, i]) for i=1:d]..., dims=[1, 2])
    @test diamond_distance(DynamicalMatrix(J1, d, d), DynamicalMatrix(J2, d, d)) ≈ 2 atol=1e-5
end

@testset "diamond distance symmetry" begin
    p=0.2
    AD=KrausOperators([[1 0; 0 sqrt(p)], [0 sqrt(1-p); 0 0]])
    AD2=KrausOperators([[1 0; 0 sqrt(2p)], [0 sqrt(1-2p); 0 0]])
    @test diamond_distance(AD,AD2) ≈ diamond_distance(AD2,AD) atol=1e-5
end

end

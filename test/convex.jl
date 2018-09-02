@testset "Convex" begin

@testset "diamond norm" begin
    srand(42)
    c = ChoiJamiolkowskiMatrices(3, 4)
    Φ = rand(c)
    @test norm_diamond(Φ) ≈ 1 atol=1e-6
end

@testset "diamond distance" begin
    d = 4
    U1 = qft(d)
    U2 = Matrix{ComplexF64}(I, d, d)
    J1 = cat([1, 2], [proj(U1[:, i]) for i=1:d]...)
    J2 = cat([1, 2], [proj(U2[:, i]) for i=1:d]...)
    @test diamond_distance(DynamicalMatrix(J1,d,d), DynamicalMatrix(J2,d,d)) ≈ 2 atol=1e-6
end

end

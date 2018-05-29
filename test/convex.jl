@testset "Convex" begin

@testset "diamond_norm" begin
    srand(42)
    J = random_dynamical_matrix(3, 4)
    @test diamond_norm(J, 3, 4) ≈ 1 atol=1e-6
end

@testset "diamond_distance" begin
    d = 4
    U1 = qft(d)
    U2 = eye(ComplexF64, d)
    J1 = cat([1, 2], [proj(U1[:, i]) for i=1:d]...)
    J2 = cat([1, 2], [proj(U2[:, i]) for i=1:d]...)
    @test diamond_distance(J1, J2, d) ≈ 2 atol=1e-6
end

end

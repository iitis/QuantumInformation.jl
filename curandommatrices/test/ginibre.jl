Random.seed!(42)

@testset "GinibreEnsemble" begin
    g = GinibreEnsemble{1}(10, 20)
    z = rand(g)
    @test eltype(z) <: Real
    @test size(z) == (10, 20)
    @test GinibreEnsemble(10, 20) == GinibreEnsemble{2}(10, 20)
    @test GinibreEnsemble(10) == GinibreEnsemble{2}(10, 10)

    @test_throws ArgumentError GinibreEnsemble{4}(11, 21)
    g = GinibreEnsemble{4}(10, 20)
    z = rand(g)
    @test size(z) == (20, 40)
    @test eltype(z) == ComplexF64

@testset "_qr_fix" begin
    a = rand(2, 2)
    u1 = RandomMatrices._qr_fix(a)
    u2 = RandomMatrices._qr_fix!(a)
    @test norm(u1 - u2) ≈ 0
end
end

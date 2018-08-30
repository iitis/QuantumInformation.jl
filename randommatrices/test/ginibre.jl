Random.seed!(42)

@testset "GinibreEnsemble" begin
    g = GinibreEnsemble{1}(10, 20)
    z = rand(g)
    @test eltype(z) <: Real
    @test size(z) == (10, 20)

    @test_throws ArgumentError GinibreEnsemble{4}(11, 21)
    g = GinibreEnsemble{4}(10, 20)
    @test_throws Exception rand(g)
end

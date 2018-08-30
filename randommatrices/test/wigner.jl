Random.seed!(42)

@testset "WignerEnsemble" begin
    w = WignerEnsemble{1}(10)
    z = rand(w)

    @test size(z) == (10, 10)
    @test eltype(z) <: Real
    @test ishermitian(z)
end

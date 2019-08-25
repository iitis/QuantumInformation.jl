Random.seed!(42)

@testset "WishartEnsemble" begin
    w = WishartEnsemble{1, 0.1}(10)
    z = rand(w)
    ev = eigvals(z)

    @test size(z) == (10, 10)
    @test eltype(z) <: Real
    @test length(ev[ev.>1e-5]) == 1
    @test all(ev[ev.>1e-5] .> 0)

    @test WishartEnsemble{1}(10) == WishartEnsemble{1, 1}(10)
    @test WishartEnsemble(10) == WishartEnsemble{2, 1}(10)
end

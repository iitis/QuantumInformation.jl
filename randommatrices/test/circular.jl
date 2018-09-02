Random.seed!(42)

@testset "CircularEnsemble" begin
    @testset "CUE" begin
        n=10
        c = CUE(n)
        u = rand(c)
        @test norm(u*u' - I) ≈ 0 atol=1e-13

        n = 100
        c = CUE(n)
        steps = 100
        r = zeros(steps, n)

        for i=1:steps
            u = rand(c)
            r[i, :] = angle.(eigvals(u))
        end
        r = vec(r)
        h = normalize(fit(Histogram, r, weights(ones(size(r))), -π:0.1π:π, closed=:left))
        @test all(isapprox.(h.weights, 1/2π, atol=0.01))
    end

    @testset "COE" begin
        n = 10
        c = COE(n)
        o = rand(c)
        @test norm(o*o' - I) ≈ 0 atol=1e-13

        n = 100
        c = COE(n)
        steps = 100
        r = zeros(steps, n)
        for i=1:steps
            o = rand(c)
            r[i, :] = angle.(eigvals(o))
        end
        r = vec(r)
        h = normalize(fit(Histogram, r, weights(ones(size(r))), -π:0.1π:π, closed=:left))
        @test all(isapprox.(h.weights, 1/2π, atol=0.1))
    end

    @testset "CircularRealEnsemble" begin
        c = CircularRealEnsemble(10)
        o = rand(c)
        @test size(o) == (10, 10)
        @test eltype(o) <: Real
    end
end

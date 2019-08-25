Random.seed!(42)

@testset "CircularEnsemble" begin
    @testset "CUE" begin
        n = 10
        c = CUE(n)
        u = CuRandomMatrices.curand(c)
        @test norm(u*u' - I) ≈ 0 atol=1e-5

        n = 100
        c = CUE(n)
        steps = 100
        r = zeros(steps, n)

        for i=1:steps
            u = collect(CuRandomMatrices.curand(c))
            r[i, :] = angle.(eigvals(u))
        end
        r = vec(r)
        h = normalize(fit(Histogram, r, weights(ones(size(r))), -π:0.1π:π, closed=:left))
        @test all(isapprox.(h.weights, 1/2π, atol=0.01))
    end

    @testset "COE" begin
        n = 10
        c = COE(n)
        o = CuRandomMatrices.curand(c)
        @test norm(o*o' - I) ≈ 0 atol=1e-5

        n = 100
        c = COE(n)
        steps = 100
        r = zeros(steps, n)
        for i=1:steps
            o = collect(CuRandomMatrices.curand(c))
            r[i, :] = angle.(eigvals(o))
        end
        r = vec(r)
        h = normalize(fit(Histogram, r, weights(ones(size(r))), -π:0.1π:π, closed=:left))
        @test all(isapprox.(h.weights, 1/2π, atol=0.1))
    end

    @testset "CircularRealEnsemble" begin
        c = CircularRealEnsemble(10)
        o = CuRandomMatrices.curand(c)
        @test size(o) == (10, 10)
        @test eltype(o) <: Real
    end
end

    @testset "HaarIsometry" begin
        idim = 2
        odim = 3
        c = HaarIsometry(idim, odim)
        u = CuRandomMatrices.curand(c)
        @test size(u) == (odim, idim)
        @test isapprox(norm(u'*u - I), 0, atol=1e-6)
        @test_throws ArgumentError HaarIsometry(odim, idim)

    @testset "CSE" begin
        n = 10
        c = CSE(n)
        o = CuRandomMatrices.curand(c)
        @test norm(o*o' - I) ≈ 0 atol=1e-5
        @test size(o) == (n, n)
    end

    @testset "Circular quaternion ensemble" begin
        c = CircularQuaternionEnsemble(10)
        u = CuRandomMatrices.curand(c)
        @test size(u) == (20, 20)
        @test isapprox(norm(u'*u - I), 0, atol=1e-5)
    end
end
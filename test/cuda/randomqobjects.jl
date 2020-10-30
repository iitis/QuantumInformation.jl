@testset "curandomqobjects" begin
    @testset "HaarKet" begin
        d = 10
        h1 = HaarKet{1}(d)
        h2 = HaarKet{2}(d)

        ϕ = curand(h1)
        ψ = curand(h2)

        @test length(ϕ) == d
        @test length(ψ) == d
        @test typeof(ϕ) == CuVector{Float32}
        @test typeof(ψ) == CuVector{ComplexF32}
        @test sum(abs2.(ϕ)) ≈ 1. atol=1e-6
        @test sum(abs2.(ψ)) ≈ 1. atol=1e-6

        h2 = HaarKet(d)
        ϕ = curand(h2)

        @test length(ϕ) == d
        @test typeof(ϕ) == CuVector{ComplexF32}
        @test sum(abs2.(ϕ)) ≈ 1. atol=1e-6
    end

    @testset "HilbertSchmidtStates" begin
        d = 10
        hs = HilbertSchmidtStates{1, 0.1}(d)
        ρ = curand(hs)

        @test size(ρ) == (d, d)
        @test typeof(ρ) == CuMatrix{Float32}
        @test norm(ρ - ρ') ≈ 0. atol=1e-13 # is close to hermitian
        @test tr(ρ) ≈ 1. atol=1e-6

        hs = HilbertSchmidtStates{2, 0.1}(d)
        ρ = curand(hs)

        @test size(ρ) == (d, d)
        @test typeof(ρ) == CuMatrix{ComplexF32}
        @test norm(ρ - ρ') ≈ 0. atol=1e-13 # is close to hermitian
        @test tr(ρ) ≈ 1. atol=1e-6

        @test HilbertSchmidtStates{1}(d) == HilbertSchmidtStates{1, 1}(d)
        @test HilbertSchmidtStates(d) == HilbertSchmidtStates{2, 1}(d)
    end

    @testset "ChoiJamiolkowskiMatrices" begin
        idim = 5
        odim = 6
        c = ChoiJamiolkowskiMatrices{1, 0.1}(idim, odim)
        j = curand(c)
        @test ptrace(j.CuMatrix, [odim, idim], [1]) ≈ I atol=1e-13
        @test tr(j.CuMatrix) ≈ idim atol=1e-13
        @test typeof(j.CuMatrix) == CuMatrix{Float32}

        c = ChoiJamiolkowskiMatrices{2, 0.1}(idim, odim)
        j = curand(c)
        @test ptrace(j.CuMatrix, [odim, idim], [1]) ≈ I atol=1e-13
        @test tr(j.CuMatrix) ≈ idim atol=1e-13
        @test typeof(j.CuMatrix) == CuMatrix{ComplexF32}

        @test ChoiJamiolkowskiMatrices{1}(idim) == ChoiJamiolkowskiMatrices{1, 1}(idim, idim)
        @test ChoiJamiolkowskiMatrices{1}(idim, odim) == ChoiJamiolkowskiMatrices{1, 1}(idim, odim)
        @test ChoiJamiolkowskiMatrices(idim, odim) == ChoiJamiolkowskiMatrices{2, 1}(idim, odim)
        @test ChoiJamiolkowskiMatrices(idim) == ChoiJamiolkowskiMatrices{2, 1}(idim, idim)
    end

    @testset "HaarPOVMs" begin
        idim = 2
        odim = 3
        c = HaarPOVM(idim, odim)
        @test_throws ArgumentError HaarPOVM(odim, idim)

        p = curand(c)
        @test norm(sum(p.matrices) - I) ≈ 0  atol=1e-8
    end

    @testset "VonNeumannPOVMs" begin
        d = 3
        c = VonNeumannPOVM(d)

        p = curand(c)
        @test norm(sum(p.matrices) - I) ≈ 0  atol=1e-8
        @test length(p.matrices) == d
    end

    @testset "WishartPOVMs" begin
        idim = 2
        odim = 3
        c = WishartPOVM(idim, odim)

        p = curand(c)
        @test norm(sum(p.matrices) - I) ≈ 0  atol=1e-8
    end
end

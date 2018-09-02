Random.seed!(42)

@testset "randomqobjects" begin
    @testset "HaarKet" begin
        d = 10
        h1 = HaarKet{1}(d)
        h2 = HaarKet{2}(d)

        ϕ = rand(h1)
        ψ = rand(h2)

        @test length(ϕ) == d
        @test length(ψ) == d
        @test typeof(ϕ) == Vector{Float64}
        @test typeof(ψ) == Vector{ComplexF64}
        @test sum(abs2.(ϕ)) ≈ 1. atol=1e-15
        @test sum(abs2.(ψ)) ≈ 1. atol=1e-15

        h2 = HaarKet(d)
        ϕ = rand(h2)

        @test length(ϕ) == d
        @test typeof(ϕ) == Vector{ComplexF64}
        @test sum(abs2.(ϕ)) ≈ 1. atol=1e-15
    end

    @testset "HilbertSchmidtStates" begin
        d = 10
        hs = HilbertSchmidtStates{1, 0.1}(d)
        ρ = rand(hs)

        @test size(ρ) == (d, d)
        @test typeof(ρ) == Matrix{Float64}
        @test norm(ρ - ρ') ≈ 0. atol=1e-13 # is close to hermitian
        @test tr(ρ) ≈ 1. atol=1e-15

        hs = HilbertSchmidtStates{2, 0.1}(d)
        ρ = rand(hs)

        @test size(ρ) == (d, d)
        @test typeof(ρ) == Matrix{ComplexF64}
        @test norm(ρ - ρ') ≈ 0. atol=1e-13 # is close to hermitian
        @test tr(ρ) ≈ 1. atol=1e-15

        @test HilbertSchmidtStates{1}(d) == HilbertSchmidtStates{1, 1}(d)
        @test HilbertSchmidtStates(d) == HilbertSchmidtStates{2, 1}(d)
    end

    @testset "ChoiJamiolkowskiMatrices" begin
        idim = 5
        odim = 6
        c = ChoiJamiolkowskiMatrices{1, 0.1}(idim, odim)
        j = rand(c)
        @test ptrace(j.matrix, [odim, idim], [1]) ≈ I atol=1e-13
        @test tr(j.matrix) ≈ idim atol=1e-13
        @test typeof(j.matrix) == Matrix{Float64}

        c = ChoiJamiolkowskiMatrices{2, 0.1}(idim, odim)
        j = rand(c)
        @test ptrace(j.matrix, [odim, idim], [1]) ≈ I atol=1e-13
        @test tr(j.matrix) ≈ idim atol=1e-13
        @test typeof(j.matrix) == Matrix{ComplexF64}

        @test ChoiJamiolkowskiMatrices{1}(idim) == ChoiJamiolkowskiMatrices{1, 1}(idim, idim)
        @test ChoiJamiolkowskiMatrices{1}(idim, odim) == ChoiJamiolkowskiMatrices{1, 1}(idim, odim)
        @test ChoiJamiolkowskiMatrices(idim, odim) == ChoiJamiolkowskiMatrices{2, 1}(idim, odim)
        @test ChoiJamiolkowskiMatrices(idim) == ChoiJamiolkowskiMatrices{2, 1}(idim, idim)
    end
end

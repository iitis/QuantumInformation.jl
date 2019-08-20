@testset "Channels" begin

include("test_channels.jl")

@testset "KrausOperators" begin
    @testset "construction" begin
        kl = [[1 0; 0 1], [1 0 0; 1 0 0; 0 0 1],[1 0; 0 1]]
        @test_throws ArgumentError KrausOperators(kl)
    end

    @testset "iscptp" begin
        for kraus_list in kraus_set
            Φ = KrausOperators(kraus_list)
            @test iscptp(Φ) == true
        end
    end

    @testset "convert to SuperOperator" begin
        ket0 = ket(1, 2)
        ket1 = ket(2, 2)
        ko = KrausOperators(kraus_list_u)
        Φ = convert(SuperOperator{Matrix{ComplexF64}}, ko)
        @test Φ(proj(ket1)) - proj(ket0) ≈ zero(proj(ket0))

        for kraus_list in kraus_set
            r, c = size(kraus_list[1])
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            Φ = convert(SuperOperator{T}, ko)
            @test ispositive(reshuffle(Φ.matrix, [r r; c c])) == true
        end
    end

    @testset "convert to Stinespring" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            Φ = convert(Stinespring{T}, ko)
            u = Φ.matrix
            @test isidentity(u'*u) == true
        end
    end

    @testset "convert to DynamicalMatrix" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            Φ = convert(DynamicalMatrix{T}, ko)
            r, c = size(kraus_list[1])
            @test isidentity(ptrace(Φ.matrix, [r, c], 1)) == true
        end
    end
end

@testset "SuperOperator" begin
    @testset "construction from function" begin
        ρ = [0.25 0.25im; -0.25im 0.75]
        t = hcat([ComplexF64[0.25, 0.25im, -0.25im, 0.75] for i=1:4]...) #stack res ρ
        m = SuperOperator{Matrix{ComplexF64}}(x -> ρ, 2, 2).matrix
        @test norm(t-m) ≈ 0. atol=1e-15
    end

    @testset "convert to KrausOperators" begin
        for kraus_list in kraus_set
            ko1 = KrausOperators(kraus_list)
            T = typeof(ko1.matrices[1])
            Φ1 = convert(SuperOperator{T}, ko1)
            ko2 = convert(KrausOperators{T}, Φ1)
            Φ2 = convert(SuperOperator{T}, ko2)
            @test iscptp(ko2)
            @test Φ1.matrix ≈ Φ2.matrix
        end
    end
    @testset "convert to DynamicalMatrix" begin
        for kraus_list in kraus_set
            r, c = size(kraus_list[1])
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            s = convert(SuperOperator{T}, ko)
            r1 = reshuffle(s.matrix, [r r; c c])
            Φ = convert(DynamicalMatrix{T}, s)
            r2 = Φ.matrix
            @test r1 ≈ r2
            @test ispositive(r2)
            @test isidentity(ptrace(r2, [r, c], 1)) == true
        end
    end
    @testset "convert to Stinespring" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            s = convert(SuperOperator{T}, ko)
            Φ = convert(Stinespring{T}, s)
            u = Φ.matrix
            @test isidentity(u'*u) == true
        end
    end
end

@testset "DynamicalMatrix" begin
    @testset "convert to KrausOperators" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            d = convert(DynamicalMatrix{T}, ko)
            kl = convert(KrausOperators{T}, d)
            @test iscptp(kl) == true
        end
    end
    @testset "convert to Stinespring" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            r = convert(DynamicalMatrix{T}, ko)
            Φ = convert(Stinespring{T}, r)
            u = Φ.matrix
            @test isidentity(u'*u) == true
        end
    end
    @testset "convert to SuperOperator" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            r1 = convert(DynamicalMatrix{T}, ko)
            s1 = convert(SuperOperator{T}, ko)
            s2 = convert(SuperOperator{T}, r1)
            @test s1.matrix ≈ s2.matrix
        end
    end
    @test_throws ArgumentError DynamicalMatrix(rand(4,5), 4, 5)
end

@testset "UnitaryChannel" begin
    @test_throws ArgumentError UnitaryChannel(ones(4, 5))
    @test_throws ArgumentError UnitaryChannel(ones(4, 4), 4, 5)

    c = UnitaryChannel(Diagonal(ComplexF64[1 -1.0im]))
    @test c isa UnitaryChannel{<:Diagonal}
end

@testset "POVMMeasurement" begin
    @testset "convert from KrausOperators" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            p = POVMMeasurement{T}(ko.matrices)
            @test !ispovm(p)
        end
    end
end

@testset "Channels applications" begin
    α = 0.25
    K₁ = ComplexF64[0 sqrt(α); 0 0; 0 0]
    K₂ = ComplexF64[1 0; 0 0; 0 sqrt(1 - α)]
    kl = Matrix{ComplexF64}[K₁, K₂]
    ρ = [0.25 0.25im; -0.25im 0.75]

    ξ = ComplexF64[0.4375 0 0.21650635im; 0 0 0; -0.21650635im 0 0.5625]

    @testset "KrausOperators" begin
        σ = KrausOperators(kl)(ρ)
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-8
    end

    @testset "DynamicalMatrix" begin
        ko = KrausOperators(kl)
        T = typeof(ko.matrices[1])
        Φ = convert(DynamicalMatrix{T}, ko)
        σ = Φ(ρ)
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-8
    end

    @testset "SuperOperator" begin
        ko = KrausOperators(kl)
        T = typeof(ko.matrices[1])
        Φ = convert(SuperOperator{T}, ko)
        σ = Φ(ρ)
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-8
    end

    @testset "Stinespring" begin
        ko = KrausOperators(kl)
        T = typeof(ko.matrices[1])
        Φ = convert(Stinespring{T}, ko)
        σ = Φ(ρ)
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-8
    end
end
end

@testset "represent" begin
    for kraus_list in kraus_set
        Φ = KrausOperators(kraus_list)
        @test represent(Φ) == kraus_list
    end

    @test represent(DynamicalMatrix(J_random, 3, 3)) == J_random
end

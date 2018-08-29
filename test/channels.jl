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
        ket0 = ket(0, 2)
        ket1 = ket(1, 2)
        @test SuperOperator{Matrix{ComplexF64}}(KrausOperators(kraus_list_u))(proj(ket1)) - proj(ket0) ≈ zero(proj(ket0))

        for kraus_list in kraus_set
            r, c = size(kraus_list[1])
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            s = SuperOperator{T}(ko)
            @test ispositive(reshuffle(s.matrix, [r r; c c])) == true
        end
    end

    @testset "convert to Stinespring" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            st = Stinespring{T}(ko)
            u = st.matrix
            @test isidentity(u'*u) == true
        end
    end

    @testset "convert to DynamicalMatrix" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            st = DynamicalMatrix{T}(ko)
            r, c = size(kraus_list[1])
            @test isidentity(ptrace(st.matrix, [r, c], 1)) == true
        end
    end
end

@testset "SuperOperator" begin
    @testset "construction from function" begin
        ρ = [0.25 0.25im; -0.25im 0.75]
        t = hcat([ComplexF64[0.25, 0.25im, -0.25im, 0.75] for i=1:4]...) #stack res ρ
        @test_throws ErrorException m = SuperOperator{Matrix{ComplexF64}}(x -> ρ, 2, 2).matrix
        @test_broken norm(t-m) ≈ 0. atol=1e-15
    end

    @testset "convert to KrausOperators" begin
        for kraus_list in kraus_set
            ko1 = KrausOperators(kraus_list)
            T = typeof(ko1.matrices[1])
            s1 = SuperOperator{T}(ko1)
            ko2 = KrausOperators{T}(s1)
            s2 = SuperOperator{T}(ko2)
            @test iscptp(ko2)
            @test s1.matrix ≈ s2.matrix
        end
    end
    @testset "convert to DynamicalMatrix" begin
        for kraus_list in kraus_set
            r, c = size(kraus_list[1])
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            s = SuperOperator{T}(ko)
            r1 = reshuffle(s.matrix, [r r; c c])
            r2 = DynamicalMatrix{T}(s).matrix
            @test r1 ≈ r2
            @test ispositive(r2)
            @test isidentity(ptrace(r2, [r, c], 1)) == true
        end
    end
    @testset "convert to Stinespring" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            s = SuperOperator{T}(ko)
            u = Stinespring{T}(s).matrix
            @test isidentity(u'*u) == true
        end
    end
end

@testset "DynamicalMatrix" begin
    @testset "convert to KrausOperators" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            r = DynamicalMatrix{T}(ko)
            kl = KrausOperators{T}(r)
            @test iscptp(kl) == true
        end
    end
    @testset "convert to Stinespring" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            r = DynamicalMatrix{T}(ko)
            u = Stinespring{T}(r).matrix
            @test isidentity(u'*u) == true
        end
    end
    @testset "convert to SuperOperator" begin
        for kraus_list in kraus_set
            ko = KrausOperators(kraus_list)
            T = typeof(ko.matrices[1])
            r1 = DynamicalMatrix{T}(ko)
            s1 = SuperOperator{T}(ko).matrix
            s2 = SuperOperator{T}(r1).matrix
            @test s1 ≈ s2
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
        σ = DynamicalMatrix{T}(ko)(ρ)
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-8
    end

    @testset "SuperOperator" begin
        ko = KrausOperators(kl)
        T = typeof(ko.matrices[1])
        σ = SuperOperator{T}(ko)(ρ)
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-8
    end

    @testset "Stinespring" begin
        ko = KrausOperators(kl)
        T = typeof(ko.matrices[1])
        σ = Stinespring{T}(ko)(ρ)
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-8
    end
end
end

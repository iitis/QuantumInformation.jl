sparsify_kraus(kraus_list) = [sparse(k) for k in kraus_list]


@testset "Channels" begin

include("test_channels.jl")

@testset "SuperOperator construction from function" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    t = hcat([ComplexF64[0.25, 0.25im, -0.25im, 0.75] for i=1:4]...) #stack res ρ
    @test_throws ErrorException m = SuperOperator{Matrix{ComplexF64}}(x -> ρ, 2, 2).matrix
    @test_broken norm(t-m) ≈ 0. atol=1e-15
end


@testset "KrausOperators" begin
    @testset "KrausOperators construction" begin
        kl = [[1 0; 0 1], [1 0 0; 1 0 0; 0 0 1],[1 0; 0 1]]
        @test_throws ArgumentError KrausOperators(kl)
    end

    @testset "KrausOperators iscptp" begin
        for kraus_list in kraus_set
            Φ = KrausOperators(kraus_list)
            @test iscptp(Φ) == true
        end
        for kraus_list in kraus_set
            Φ = KrausOperators(sparsify_kraus(kraus_list))
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

if false
@testset "superoperator" begin
    @testset "superoperator_to_kraus" begin

    for kraus_list in kraus_set
        s1 = kraus_to_superoperator(kraus_list)
        kl = superoperator_to_kraus(s1)
        s2 = kraus_to_superoperator(kl)
        @test s1 ≈ s2
    end
    for kraus_list in kraus_set
        s1 = kraus_to_superoperator(kraus_list)
        kl = superoperator_to_kraus(sparse(s1))
        s2 = kraus_to_superoperator(kl)
        @test s1 ≈ s2
    end

    end
    @testset "superoperator_to_dynamical_matrix" begin
    for kraus_list in kraus_set
        r, c = size(kraus_list[1])
        s = kraus_to_superoperator(kraus_list)
        r1 = reshuffle(s, [r r; c c])
        r2 = superoperator_to_dynamical_matrix(s)
        @test r1 ≈ r2
        @test ispositive(r2)
        @test isidentity(ptrace(r2, [r, c], 1)) == true
    end

    for kraus_list in kraus_set
        r, c = size(kraus_list[1])
        s = kraus_to_superoperator(sparsify_kraus(kraus_list))
        r1 = reshuffle(s, [r r; c c])
        r2 = superoperator_to_dynamical_matrix(s)
        @test r1 ≈ r2
        @test issparse(r2)
        @test ispositive(r2)
        @test isidentity(ptrace(r2, [r, c], 1)) == true
    end


    end
    @testset "superoperator_to_stinespring" begin
        for kraus_list in kraus_set
            s = kraus_to_superoperator(kraus_list)
            u = superoperator_to_stinespring(s)
            @test isidentity(u'*u) == true
        end
        for kraus_list in kraus_set
            s = kraus_to_superoperator(kraus_list)
            u = superoperator_to_stinespring(sparse(s))
            @test isidentity(u'*u) == true
        end
    end
end

@testset "dynamical matrix" begin
    @testset "dynamical_matrix_to_kraus" begin
        for kraus_list in kraus_set
            r = kraus_to_dynamical_matrix(kraus_list)
            rows, cols = size(kraus_list[1])
            kl = dynamical_matrix_to_kraus(r, rows, cols)
            @test kraus_is_CPTP(kl) == true
        end
        for kraus_list in kraus_set
            r = kraus_to_dynamical_matrix(kraus_list)
            rows, cols = size(kraus_list[1])
            kl = dynamical_matrix_to_kraus(sparse(r), rows, cols)
            @test kraus_is_CPTP(kl) == true
        end
    end
    @testset "dynamical_matrix_to_stinespring" begin
        for kraus_list in kraus_set
            rows, cols = size(kraus_list[1])
            r = kraus_to_dynamical_matrix(kraus_list)
            u = dynamical_matrix_to_stinespring(r, rows, cols)
            @test isidentity(u'*u) == true
        end
        for kraus_list in kraus_set
            r, c = size(kraus_list[1])
            rows, cols = size(kraus_list[1])
            r = kraus_to_dynamical_matrix(kraus_list)
            u = dynamical_matrix_to_stinespring(sparse(r), rows, cols)
            @test isidentity(u'*u) == true
        end
    end
    @testset "dynamical_matrix_to_superoperator" begin
        for kraus_list in kraus_set
            r1 = kraus_to_dynamical_matrix(kraus_list)
            s = kraus_to_superoperator(kraus_list)
            rows, cols = size(kraus_list[1])
            r2 = dynamical_matrix_to_superoperator(s, rows, cols)
            @test r1 ≈ r2
        end
        for kraus_list in kraus_set
            r1 = kraus_to_dynamical_matrix(kraus_list)
            s = kraus_to_superoperator(kraus_list)
            rows, cols = size(kraus_list[1])
            r2 = dynamical_matrix_to_superoperator(sparse(s), rows, cols)
            @test r1 ≈ r2
        end
    end
end

@testset "apply_channel" begin
    α = 0.25
    K₁ = ComplexF64[0 sqrt(α); 0 0]
    K₂ = ComplexF64[1 0; 0 sqrt(1 - α)]
    kl = Matrix{ComplexF64}[K₁, K₂]
    ρ = [0.25 0.25im; -0.25im 0.75]
    ξ = ComplexF64[1/4 + 3/16 sqrt(3/4)*1im/4; -sqrt(3/4)*1im/4 9/16]

    @testset "apply_channel_kraus" begin
        σ = apply_channel_kraus(kl, ρ)
        @test trace(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-15
    end

    @testset "apply_channel_dynamical_matrix" begin
        r = kraus_to_dynamical_matrix(kl)
        σ = apply_channel_dynamical_matrix(r, ρ)
        @test trace(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-15
    end

    @testset "apply_channel_superoperator" begin
        s = kraus_to_superoperator(kl)
        σ = apply_channel_superoperator(s, ρ)
        @test trace(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-15
    end

    @testset "apply_channel_stinespring" begin
        u = kraus_to_stinespring(kl)
        dims = [size(kl[1])...]
        σ = apply_channel_stinespring(u, ρ, dims)
        @test trace(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-15
    end
end
end # if false
end

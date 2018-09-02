sparsify_kraus(kraus_list) = [sparse(k) for k in kraus_list]


@testset "Channels" begin

include("test_channels.jl")

@testset "channel_to_superoperator" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    T = hcat([ComplexF64[0.25, 0.25im, -0.25im, 0.75] for i=1:4]...) #stack res ρ
    M = channel_to_superoperator(x -> ρ, 2)
    @test norm(T-M) ≈ 0. atol=1e-15
end

@testset "kraus_check_size" begin
    for kraus_list in kraus_set
        @test kraus_check_size(kraus_list) == nothing
    end
    for kraus_list in kraus_set
        @test kraus_check_size(sparsify_kraus(kraus_list)) == nothing
    end
    kl = [[1 0; 0 1], [1 0 0; 1 0 0; 0 0 1],[1 0; 0 1]]
    @test_throws ArgumentError kraus_check_size(kl)
end

@testset "kraus" begin
    @testset "kraus_is_CPTP" begin
        for kraus_list in kraus_set
            @test kraus_is_CPTP(kraus_list) == true
        end
        for kraus_list in kraus_set
            @test kraus_is_CPTP(sparsify_kraus(kraus_list)) == true
        end
    end

    @testset "kraus_to_superoperator" begin
        ket0 = ket(0, 2)
        ket1 = ket(1, 2)
        @test kraus_to_superoperator(kraus_list_u)*vec(proj(ket1)) ≈ vec(proj(ket0))

        ket0 = ket(SparseVector{ComplexF64}, 0, 2)
        ket1 = ket(SparseVector{ComplexF64}, 1, 2)
        @test kraus_to_superoperator(sparsify_kraus(kraus_list_u))*vec(proj(ket0)) ≈ vec(proj(ket1))

        for kraus_list in kraus_set
            r, c = size(kraus_list[1])
            @test ispositive(reshuffle(kraus_to_superoperator(kraus_list), [r r; c c])) == true
        end
        for kraus_list in kraus_set
            r, c = size(kraus_list[1])
            @test ispositive(reshuffle(kraus_to_superoperator(sparsify_kraus(kraus_list)), [r r; c c])) == true
        end
    end

    @testset "kraus_to_stinespring" begin
        for kraus_list in kraus_set
            u = kraus_to_stinespring(kraus_list)
            @test isidentity(u'*u)==true
        end
        for kraus_list in kraus_set
            u = kraus_to_stinespring(sparsify_kraus(kraus_list))
            @test isidentity(u'*u)==true
        end
    end

    @testset "kraus_to_dynamical_matrix" begin
        for kraus_list in kraus_set
            h = kraus_to_dynamical_matrix(kraus_list)
            r, c = size(kraus_list[1])
            @test isidentity(ptrace(h, [r, c], 1))==true
        end
        for kraus_list in kraus_set
            h = kraus_to_dynamical_matrix(sparsify_kraus(kraus_list))
            r, c = size(kraus_list[1])
            @test isidentity(ptrace(h, [r, c], 1))==true
        end
    end
end

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
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-15
    end

    @testset "apply_channel_dynamical_matrix" begin
        r = kraus_to_dynamical_matrix(kl)
        σ = apply_channel_dynamical_matrix(r, ρ)
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-15
    end

    @testset "apply_channel_superoperator" begin
        s = kraus_to_superoperator(kl)
        σ = apply_channel_superoperator(s, ρ)
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-15
    end

    @testset "apply_channel_stinespring" begin
        u = kraus_to_stinespring(kl)
        dims = [size(kl[1])...]
        σ = apply_channel_stinespring(u, ρ, dims)
        @test tr(σ) ≈ 1. atol=1e-15
        @test ishermitian(σ)
        @test σ ≈ ξ atol=1e-15
    end
end

end

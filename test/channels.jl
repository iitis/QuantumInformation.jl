@testset "Channels" begin

@testset "kraus" begin

    kraus_list_u = [sx]

    α = 0.2
    kraus_list_dep = [sqrt(1-α)*eye(2), sqrt(α/3)*sx, sqrt(α/3)*sz, sqrt(α/3)*sy]

    α = 0.3
    k1 = [1 0; 0 sqrt(1-α)]
    k2 = [0 sqrt(α); 0 0]
    kraus_list_ad = [k1, k2]

    k1 = [1 0 0; 0 1 0]
    k2 = [0 0 1; 0 0 0]
    kraus_list_dim = [k1, k2]

    sparsify_kraus(kraus_list) = [sparse(k) for k in kraus_list]

    @testset "kraus_is_CPTP" begin
        @test kraus_is_CPTP(kraus_list_u) == true
        @test kraus_is_CPTP(kraus_list_dep) == true
        @test kraus_is_CPTP(kraus_list_ad) == true
        @test kraus_is_CPTP(kraus_list_dim) == true

        @test kraus_is_CPTP(sparsify_kraus(kraus_list_u)) == true
        @test kraus_is_CPTP(sparsify_kraus(kraus_list_dep)) == true
        @test kraus_is_CPTP(sparsify_kraus(kraus_list_ad)) == true
        @test kraus_is_CPTP(sparsify_kraus(kraus_list_dim)) == true
    end

    @testset "kraus_to_superoperator" begin
        ket0 = ket(0, 2)
        ket1 = ket(1, 2)
        @test kraus_to_superoperator(kraus_list_u)*vec(proj(ket1)) ≈ vec(proj(ket0))

        ket0 = ket(SparseVector{ComplexF64}, 0, 2)
        ket1 = ket(SparseVector{ComplexF64}, 1, 2)
        @test kraus_to_superoperator(sparsify_kraus(kraus_list_u))*vec(proj(ket0)) ≈ vec(proj(ket1))
    end
end
end

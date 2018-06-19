@testset "Channels" begin

include("test_channels.jl")

@testset "kraus" begin



    sparsify_kraus(kraus_list) = [sparse(k) for k in kraus_list]

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
    end
    
end
end

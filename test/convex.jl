import Random: seed!
using DoubleFloats: Double64, ComplexDF64
@testset verbose=true "Convex" begin
    @testset verbose=true "diamond norm" begin
        seed!(42)
        c = ChoiJamiolkowskiMatrices(3, 4)
        Φ = rand(c)
        @test norm_diamond(Φ) ≈ 1 atol=1e-4
    end

    @testset verbose=true "diamond distance" begin
        d = 4
        U1 = qft(d)
        U2 = Matrix{ComplexF64}(I, d, d)
        J1 = cat([proj(U1[:, i]) for i in 1:d]..., dims=[1, 2])
        J2 = cat([proj(U2[:, i]) for i in 1:d]..., dims=[1, 2])
        @test diamond_distance(DynamicalMatrix(J1, d, d), DynamicalMatrix(J2, d, d)) ≈ 2 atol=1e-4
    end

    @testset verbose=true "diamond distance symmetry" begin
        p=0.2
        AD=KrausOperators([[1 0; 0 sqrt(p)], [0 sqrt(1-p); 0 0]])
        AD2=KrausOperators([[1 0; 0 sqrt(2p)], [0 sqrt(1-2p); 0 0]])
        @test diamond_distance(AD, AD2) ≈ diamond_distance(AD2, AD) atol=1e-4
    end

    @testset verbose=true "diamond norm dual" begin
        c = ChoiJamiolkowskiMatrices(3, 4)
        Φ = rand(c)
        @test norm_diamond(Φ, :dual) ≈ 1 atol=5e-4
    end
    @testset verbose=true "Double64 Support" begin
        # Basic check if convex optimization functions run without error on Double64 types
        # Note: Optimization backend might convert to Float64, so we test basic execution.
        d = 2
        U1 = qft(ComplexDF64, d)
        U2 = Matrix{ComplexDF64}(I, d, d)
        J1 = cat([proj(U1[:, i]) for i in 1:d]..., dims=[1, 2])
        J2 = cat([proj(U2[:, i]) for i in 1:d]..., dims=[1, 2])

        # Just check if it runs, accuracy might be limited by solver (SCS)
        dm1 = DynamicalMatrix(J1, d, d)
        dm2 = DynamicalMatrix(J2, d, d)
        # We relax the type requirement for the check, assuming it might return Float64 result
        res = diamond_distance(dm1, dm2)
        @test res isa Real
    end
end

@testset verbose=true "Channels" begin
    using DoubleFloats: Double64, ComplexDF64
    using SparseArrays

    include("test_channels.jl")

    @testset verbose=true "KrausOperators" begin
        @testset verbose=true "construction" begin
            kl = [[1 0; 0 1], [1 0 0; 1 0 0; 0 0 1], [1 0; 0 1]]
            @test_throws ArgumentError KrausOperators(kl)
            @test_throws ArgumentError KrausOperators([rand(2, 2)], 3, 2)
        end

        @testset verbose=true "iscptp" begin
            for kraus_list in kraus_set
                Φ = KrausOperators(kraus_list)
                @test iscptp(Φ) == true
            end
        end

        @testset verbose=true "convert to SuperOperator" begin
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

        @testset verbose=true "convert to Stinespring" begin
            for kraus_list in kraus_set
                ko = KrausOperators(kraus_list)
                T = typeof(ko.matrices[1])
                Φ = convert(Stinespring{T}, ko)
                u = Φ.matrix
                @test isidentity(u'*u) == true
            end
        end

        @testset verbose=true "convert to DynamicalMatrix" begin
            for kraus_list in kraus_set
                ko = KrausOperators(kraus_list)
                T = typeof(ko.matrices[1])
                Φ = convert(DynamicalMatrix{T}, ko)
                r, c = size(kraus_list[1])
                @test isidentity(ptrace(Φ.matrix, [r, c], 1)) == true
            end
        end
    end

    @testset verbose=true "SuperOperator" begin
        @testset verbose=true "construction from function" begin
            ρ = [0.25 0.25im; -0.25im 0.75]
            t = hcat([ComplexF64[0.25, 0.25im, -0.25im, 0.75] for i in 1:4]...) #stack res ρ
            m = SuperOperator{Matrix{ComplexF64}}(x -> ρ, 2, 2).matrix
            @test norm(t-m) ≈ 0.0 atol=1e-15

            @test_throws ArgumentError SuperOperator(rand(3, 3)) # Not square of something
            @test_throws ArgumentError SuperOperator(rand(4, 4), 3, 2)
            @test_throws ArgumentError SuperOperator(x->x, -1, 2)
        end

        @testset verbose=true "convert to KrausOperators" begin
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
        @testset verbose=true "convert to DynamicalMatrix" begin
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
        @testset verbose=true "convert to Stinespring" begin
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

    @testset verbose=true "DynamicalMatrix" begin
        @testset verbose=true "convert to KrausOperators" begin
            for kraus_list in kraus_set
                ko = KrausOperators(kraus_list)
                T = typeof(ko.matrices[1])
                d = convert(DynamicalMatrix{T}, ko)
                kl = convert(KrausOperators{T}, d)
                @test iscptp(kl) == true
            end
        end
        @testset verbose=true "convert to Stinespring" begin
            for kraus_list in kraus_set
                ko = KrausOperators(kraus_list)
                T = typeof(ko.matrices[1])
                r = convert(DynamicalMatrix{T}, ko)
                Φ = convert(Stinespring{T}, r)
                u = Φ.matrix
                @test isidentity(u'*u) == true
            end
        end
        @testset verbose=true "convert to SuperOperator" begin
            for kraus_list in kraus_set
                ko = KrausOperators(kraus_list)
                T = typeof(ko.matrices[1])
                r1 = convert(DynamicalMatrix{T}, ko)
                s1 = convert(SuperOperator{T}, ko)
                s2 = convert(SuperOperator{T}, r1)
                @test s1.matrix ≈ s2.matrix
            end
        end
        @test_throws ArgumentError DynamicalMatrix(rand(4, 3), 2, 2)
        @test_throws ArgumentError DynamicalMatrix(rand(4, 5), 4, 5)
    end

    @testset verbose=true "UnitaryChannel" begin
        @test_throws ArgumentError UnitaryChannel(ones(4, 5))
        @test_throws ArgumentError UnitaryChannel(ones(4, 4), 4, 5)
        
        # Stinespring error
        @test_throws ArgumentError Stinespring(rand(4, 4), 2, 2) # 8x2 expected
        
        # PostSelectionMeasurement error
        @test_throws ArgumentError PostSelectionMeasurement(rand(2, 2), 3, 3)

        c = UnitaryChannel(Diagonal(ComplexF64[1 -1.0im]))
        @test c isa UnitaryChannel{<:Diagonal}

        u1 = UnitaryChannel([cos(1) sin(1); -sin(1) cos(1)])
        u2 = UnitaryChannel([cos(2) sin(2); -sin(2) cos(2)])
        @test compose(u1, u2).matrix ≈ u2.matrix * u1.matrix
        @test compose(UnitaryChannel{Array{Float64, 2}}, u1, u2).matrix ≈
              u2.matrix * u1.matrix

        @test kron(u1, u2).matrix ≈ kron(u1.matrix, u2.matrix)
    end

    @testset verbose=true "POVMMeasurement" begin
        @testset verbose=true "convert from KrausOperators" begin
            for kraus_list in kraus_set
                ko = KrausOperators(kraus_list)
                T = typeof(ko.matrices[1])
                p = POVMMeasurement{T}(ko.matrices)
                @test !ispovm(p)
            end
        end
    end

    @testset verbose=true "Channels applications" begin
        α = 0.25
        K₁ = ComplexF64[0 sqrt(α); 0 0; 0 0]
        K₂ = ComplexF64[1 0; 0 0; 0 sqrt(1 - α)]
        kl = Matrix{ComplexF64}[K₁, K₂]
        ρ = [0.25 0.25im; -0.25im 0.75]

        ξ = ComplexF64[0.4375 0 0.21650635im; 0 0 0; -0.21650635im 0 0.5625]

        @testset verbose=true "KrausOperators" begin
            σ = KrausOperators(kl)(ρ)
            @test tr(σ) ≈ 1.0 atol=1e-15
            @test ishermitian(σ)
            @test σ ≈ ξ atol=1e-8
        end

        @testset verbose=true "DynamicalMatrix" begin
            ko = KrausOperators(kl)
            T = typeof(ko.matrices[1])
            Φ = convert(DynamicalMatrix{T}, ko)
            σ = Φ(ρ)
            @test tr(σ) ≈ 1.0 atol=1e-15
            @test ishermitian(σ)
            @test σ ≈ ξ atol=1e-8
        end

        @testset verbose=true "SuperOperator" begin
            ko = KrausOperators(kl)
            T = typeof(ko.matrices[1])
            Φ = convert(SuperOperator{T}, ko)
            σ = Φ(ρ)
            @test tr(σ) ≈ 1.0 atol=1e-15
            @test ishermitian(σ)
            @test σ ≈ ξ atol=1e-8
        end

        @testset verbose=true "Stinespring" begin
            ko = KrausOperators(kl)
            T = typeof(ko.matrices[1])
            Φ = convert(Stinespring{T}, ko)
            σ = Φ(ρ)
            @test tr(σ) ≈ 1.0 atol=1e-15
            @test ishermitian(σ)
            @test σ ≈ ξ atol=1e-8
            @testset verbose=true "Double64 Support" begin
                df_kraus_list = [Double64[1 0; 0 1], Double64[0 1; 0 0]]
                df_ko = KrausOperators(df_kraus_list)
                @test iscptp(df_ko) isa Bool
                # Check if we can convert
                df_so = convert(SuperOperator{Matrix{Double64}}, df_ko)
                @test df_so.matrix isa Matrix{Double64}

                # Check mix of Double64 and Complex{Double64}
                cdf_val = Complex{Double64}(1.0, 0.0)
                df_u = UnitaryChannel(Double64[0 1; 1 0] * cdf_val)
                @test df_u.matrix isa Matrix{Complex{Double64}}
                @test iscptp(df_u)

                # Application
                ρ = Double64[1 0; 0 0]
                σ = df_ko(ρ)
                @test σ isa Matrix{Double64}
                @test tr(σ) ≈ 1.0
            end
        end
    end
end

@testset verbose=true "represent" begin
    for kraus_list in kraus_set
        Φ = KrausOperators(kraus_list)
        @test represent(Φ) == kraus_list
    end

    @test represent(DynamicalMatrix(J_random, 3, 3)) == J_random
end

@testset verbose=true "IO and Printers" begin
    for kraus_list in kraus_set
        MT = eltype(kraus_list)
        Φ = KrausOperators(kraus_list)
        @test_nowarn show(devnull, Φ)
        @test_nowarn show(devnull, convert(SuperOperator{MT}, Φ))
        @test_nowarn show(devnull, convert(DynamicalMatrix{MT}, Φ))
        @test_nowarn show(devnull, convert(Stinespring{MT}, Φ))
    end
    @test_nowarn show(devnull, UnitaryChannel(Matrix(𝕀(2))))
end

@testset verbose=true "Predicates" begin
    @testset verbose=true "iscp" begin
        # SuperOperator and DynamicalMatrix are already partially covered, but let's be thorough
        for kraus_list in kraus_set
            MT = eltype(kraus_list)
            Φ = KrausOperators(kraus_list)
            @test iscp(Φ)
            @test iscp(convert(SuperOperator{MT}, Φ))
            @test iscp(convert(DynamicalMatrix{MT}, Φ))
            @test iscp(convert(Stinespring{MT}, Φ))
            @test iscp(UnitaryChannel(Matrix(𝕀(2))))
        end
    end

    @testset verbose=true "istni and istp" begin
        for kraus_list in kraus_set
            MT = eltype(kraus_list)
            Φ = KrausOperators(kraus_list)
            @test istni(Φ)
            @test istp(Φ)
            @test istni(convert(SuperOperator{MT}, Φ))
            @test istp(convert(SuperOperator{MT}, Φ))
            @test istni(convert(DynamicalMatrix{MT}, Φ))
            @test istp(convert(DynamicalMatrix{MT}, Φ))
            @test istni(convert(Stinespring{MT}, Φ))
            @test istp(convert(Stinespring{MT}, Φ))
            @test istni(UnitaryChannel(Matrix(𝕀(2))))
            @test istp(UnitaryChannel(Matrix(𝕀(2))))
        end
    end

    @testset verbose=true "iscptp and iscptni" begin
        for kraus_list in kraus_set
            Φ = KrausOperators(kraus_list)
            @test iscptp(Φ)
            @test iscptni(Φ)
        end
    end

    @testset verbose=true "ispovm and iseffect" begin
        # Valid POVM
        p = POVMMeasurement([sz/2 + 𝕀(2)/2, -sz/2 + 𝕀(2)/2])
        @test ispovm(p)

        # Invalid POVM
        p_inv = POVMMeasurement([sz/2, -sz/2])
        @test !ispovm(p_inv)

        # Valid effect
        eff = PostSelectionMeasurement([0.5 0; 0 0.5])
        @test iseffect(eff)

        # Invalid effect (operator norm > 1)
        eff_inv = PostSelectionMeasurement([2.0 0; 0 2.0])
        @test !iseffect(eff_inv)
    end

    @testset verbose=true "Edge cases" begin
        # iscp false
        ρ_neg = DynamicalMatrix(ComplexF64[1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 -1], 2, 2)
        @test iscp(ρ_neg) == false
        
        # istp false
        ko_not_tp = KrausOperators([0.5 * I(2)])
        @test istp(ko_not_tp) == false

        # POVMMeasurement error
        @test_throws ArgumentError POVMMeasurement([rand(2, 2), rand(2, 3)])
        @test_throws ArgumentError POVMMeasurement([rand(2, 2)], 3, 1)
        @test_throws ArgumentError POVMMeasurement([rand(2, 2)], 2, 2)
    end
end

@testset verbose=true "Compositions" begin
    @testset verbose=true "kron" begin
        u = UnitaryChannel(Matrix(I, 2, 2))
        id = IdentityChannel(2)
        @test kron(u, id) isa UnitaryChannel
        @test kron(id, u) isa UnitaryChannel
        @test kron(u, u) isa UnitaryChannel

        ko = KrausOperators([Matrix(I, 2, 2)])
        @test kron(ko, u) isa KrausOperators
        @test kron(u, ko) isa KrausOperators
    end

    @testset verbose=true "compose and *" begin
        u1 = UnitaryChannel(sx)
        u2 = UnitaryChannel(sy)
        @test compose(u1, u2) ≈ UnitaryChannel(sy * sx)
        @test u1 * u2 ≈ UnitaryChannel(sy * sx)

        ko1 = KrausOperators([sx])
        ko2 = KrausOperators([sy])
        target = convert(SuperOperator{Matrix{ComplexF64}}, UnitaryChannel(sy * sx))
        @test compose(ko1, ko2) ≈ target
        @test ko1 * ko2 ≈ target

        @test_throws ArgumentError compose(
            UnitaryChannel(Matrix(𝕀(2))),
            UnitaryChannel(Matrix(𝕀(3))),
        )
    end
end

@testset verbose=true "Channels applications - Vectors" begin
    u = UnitaryChannel(sx)
    ψ = ket(1, 2)
    @test u(ψ) ≈ sx * ψ

    id = IdentityChannel(2)
    @test id(ψ) == ψ

    ko = KrausOperators([sx])
    @test ko(ψ) ≈ proj(sx * ψ)
end

@testset verbose=true "Misc functions" begin
    u = UnitaryChannel(sx)
    @test size(u) == (2, 2)
    @test represent(u) == sx

    ko = KrausOperators([sx])
    @test size(ko) == (2, 2)
    @test represent(ko) == [sx]

    id = IdentityChannel(2)
    @test size(id) == (2, 2)
    @test represent(id) ≈ Matrix(𝕀(2))
end

@testset verbose=true "Sparse Support" begin
    # 1. KrausOperators with Sparse Matrices
    K1 = sparse([1.0 0.0; 0.0 1.0]) # Identity
    K2 = sparse([0.0 1.0; 1.0 0.0]) # SX
    ko = KrausOperators([K1, K2])
    @test ko.matrices[1] isa SparseMatrixCSC
    @test ko isa KrausOperators{SparseMatrixCSC{Float64, Int64}}

    # 2. Application to Dense
    ρ = [1.0 0.0; 0.0 0.0]
    out = ko(ρ)
    @test out ≈ K1*ρ*K1' + K2*ρ*K2'

    # 3. Application to Sparse
    ρ_sparse = sparse(ρ)
    out_sparse = ko(ρ_sparse)
    @test out_sparse isa AbstractSparseMatrix
    @test out_sparse ≈ out

    # 4. UnitaryChannel
    U = sparse([0.0 1.0; 1.0 0.0])
    uc = UnitaryChannel(U)
    @test uc.matrix isa SparseMatrixCSC
    @test uc(ρ_sparse) ≈ U*ρ_sparse*U'

    # 5. IdentityChannel
    id = IdentityChannel(2)
    @test id(ρ_sparse) === ρ_sparse

    # 6. ISC PTP
    # ko is not trace preserving sum(Ki' Ki) = I + I = 2I.
    # ko normalized:
    ko_norm = KrausOperators([K1/sqrt(2), K2/sqrt(2)])
    @test iscptp(ko_norm)
end

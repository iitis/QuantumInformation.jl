@testset verbose=true "Permute systems" begin
    using DoubleFloats: Double64, ComplexDF64

    @testset verbose=true "Diagonal matrix" begin
        initial = Matrix(Diagonal([0, 0, 1, 0, 0, 1, 0, 0, 1]))
        permuted = Matrix(Diagonal([0, 0, 0, 0, 0, 0, 1, 1, 1]))
        @test sum(permutesystems(initial, [3, 3], [2, 1]) - permuted) ≈ 0.0 atol=1e-15
        @test sum(permutesystems(initial, [3, 3], [1, 2]) - initial) ≈ 0.0 atol=1e-15

        @test_throws ArgumentError permutesystems(ones(2, 3), [2, 2], [1, 2]) # Non-square
        @test_throws ArgumentError permutesystems(ones(4, 4), [2, 3], [1, 2]) # Product mismatch
        @test_throws ArgumentError permutesystems(ones(4, 4), [2, 2], [1, 3]) # Index out of range
    end

    @testset verbose=true "More complex diagonal matrix" begin
        # the following example has been generated
        # using original python implementation of permute_systems
        initial_diagonal = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        permuted_diagonal = [1, 5, 9, 13, 3, 7, 11, 15, 2, 6, 10, 14, 4, 8, 12, 16]
        initial = Matrix(Diagonal(initial_diagonal))
        initial[1, 8] = 42
        #         println(findfirst(initial, 42))
        permuted = Matrix(Diagonal(permuted_diagonal))
        permuted[1, 14] = 42
        @test sum(abs.(permutesystems(initial, [2, 2, 2, 2], [4, 3, 1, 2]) - permuted)) ≈
              0.0 atol=1e-15
    end

    @testset verbose=true "Different dimensions" begin
        # the following example has been generated
        # using original python implementation of permute_systems

        initial = diagm(0=>1:24)
        permuted = diagm(
            0=>[
                1,
                13,
                2,
                14,
                3,
                15,
                4,
                16,
                5,
                17,
                6,
                18,
                7,
                19,
                8,
                20,
                9,
                21,
                10,
                22,
                11,
                23,
                12,
                24,
            ],
        )

        @test sum(abs.(permutesystems(initial, [2, 3, 4], [2, 3, 1]) - permuted)) ≈ 0.0 atol=1e-15
    end
    @testset verbose=true "Double64 Support" begin
        # Diagonal matrix with Double64
        initial = Matrix(Diagonal(Double64[0, 0, 1, 0, 0, 1, 0, 0, 1]))
        permuted = Matrix(Diagonal(Double64[0, 0, 0, 0, 0, 0, 1, 1, 1]))
        res = permutesystems(initial, [3, 3], [2, 1])
        @test eltype(res) == Double64
        @test sum(res - permuted) ≈ 0.0 atol=1e-25
    end

    @testset verbose=true "Sparse Support" begin
        using SparseArrays
        # 1. Identity permutation
        s_id = sparse(I, 4, 4)
        @test permutesystems(s_id, [2, 2], [1, 2]) ≈ s_id
        @test permutesystems(s_id, [2, 2], [1, 2]) isa AbstractSparseMatrix

        # 2. Swap systems
        initial = sparse(Diagonal([0, 0, 1, 0, 0, 1, 0, 0, 1]))
        permuted = sparse(Diagonal([0, 0, 0, 0, 0, 0, 1, 1, 1]))
        s_out = permutesystems(initial, [3, 3], [2, 1])
        @test s_out isa AbstractSparseMatrix
        @test s_out ≈ permuted

        # 3. Different dimensions
        # From dense test:
        initial_dense = diagm(0=>1:24)
        permuted_dense = diagm(
            0=>[
                1,
                13,
                2,
                14,
                3,
                15,
                4,
                16,
                5,
                17,
                6,
                18,
                7,
                19,
                8,
                20,
                9,
                21,
                10,
                22,
                11,
                23,
                12,
                24,
            ],
        )

        s_init = sparse(initial_dense)
        s_perm = permutesystems(s_init, [2, 3, 4], [2, 3, 1])
        @test s_perm isa AbstractSparseMatrix
        @test s_perm ≈ permuted_dense

        # 4. Error checks
        @test_throws ArgumentError permutesystems(sparse(ones(2, 3)), [2, 2], [1, 2]) # Non-square
        @test_throws ArgumentError permutesystems(sparse(ones(4, 4)), [2, 3], [1, 2]) # Product mismatch
        @test_throws ArgumentError permutesystems(sparse(ones(4, 4)), [2, 2], [1, 3]) # Index out of range
    end
end

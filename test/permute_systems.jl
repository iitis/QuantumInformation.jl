@testset "Permute systems" begin

    @testset "Diagonal matrix" begin
         initial = Matrix(Diagonal([0,0,1, 0,0,1, 0,0,1]))
         permuted = Matrix(Diagonal([0,0,0,0,0,0,1,1,1]))
         @test sum(permutesystems(initial, [3,3], [2,1]) - permuted) ≈ 0. atol=1e-15
         @test sum(permutesystems(initial, [3,3], [1, 2]) - initial) ≈ 0. atol=1e-15

         @test_throws ArgumentError permutesystems(ones(2, 3), [1, 2], [1])
         @test_throws ArgumentError permutesystems(ones(2, 2), [3, 4], [2])
         @test_throws ArgumentError permutesystems(ones(4, 4), [2, 2], [3])
    end

    @testset "More complex diagonal matrix" begin
        # the following example has been generated
        # using original python implementation of permute_systems
        initial_diagonal = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
        permuted_diagonal = [1, 5, 9, 13, 3, 7, 11, 15, 2, 6, 10, 14, 4, 8, 12, 16]
        initial = Matrix(Diagonal(initial_diagonal))
        initial[1,8] = 42
#         println(findfirst(initial, 42))
        permuted = Matrix(Diagonal(permuted_diagonal))
        permuted[1, 14] = 42
        @test sum(abs.(permutesystems(initial, [2,2,2,2], [4,3,1,2]) - permuted)) ≈ 0. atol=1e-15

    end
end

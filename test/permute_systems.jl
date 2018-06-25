@testset "Permute systems" begin

    @testset "Diagonal matrix" begin
         initial = diagm([0,0,1, 0,0,1, 0,0,1])
         permuted = diagm([0,0,0,0,0,0,1,1,1])
         @test sum(permutesystems(initial, [3,3], [2,1]) - permuted) ≈ 0. atol=1e-15
         @test sum(permutesystems(initial, [3,3], [1, 2]) - initial) ≈ 0. atol=1e-15

         @test_throws ArgumentError permutesystems(ones(2, 3), [1, 2], [1])
         @test_throws ArgumentError permutesystems(ones(2, 2), [3, 4], [2])
         @test_throws ArgumentError permutesystems(ones(4, 4), [2, 2], [3])
    end
end

@testset verbose=true "Partial transpose" begin
using DoubleFloats: Double64, ComplexDF64

@testset verbose=true "Dense matrices" begin
  ρ =  [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]
  trans1 = [1 2 9 10; 5 6 13 14; 3 4 11 12; 7 8 15 16]
  trans2 = [1 5 3 7; 2 6 4 8; 9 13 11 15; 10 14 12 16]
  
  res1 = ptranspose(ρ, [2, 2], [1])
  res2 = ptranspose(ρ, [2, 2], [2])

  @test norm(res1 - trans1) ≈ 0. atol=1e-15
  @test norm(res2 - trans2) ≈ 0. atol=1e-15

  @test_throws ArgumentError ptranspose(ones(2, 3), [2, 2], 1)
  @test_throws ArgumentError ptranspose(ones(4, 4), [2, 3], 1)
  @test_throws ArgumentError ptranspose(ones(4, 4), [2, 2], 3)

  # Added for coverage
  @testset verbose=true "Varied dimensions" begin
      ρ3 = rand(8, 8)
      # 2x2x2 system
      res3 = ptranspose(ρ3, [2, 2, 2], [1, 3])
      @test size(res3) == (8, 8)
      
      # 2x4 system
      ρ4 = rand(8, 8)
      @test size(ptranspose(ρ4, [2, 4], [1])) == (8, 8)
      @test size(ptranspose(ρ4, [2, 4], [2])) == (8, 8)
  end

  @testset verbose=true "Internal _ptranspose" begin
      ρ = Matrix{ComplexF64}(I, 4, 4)
      # Test the internal helper which uses @cast
      res = QuantumInformation._ptranspose(ρ, [2, 2], [1])
      @test size(res) == (4, 4)
      @test res ≈ ptranspose(ρ, [2, 2], [1])
  end
end

@testset verbose=true "Double64 Support" begin
    # Use same matrix as dense test but cast to Double64
    ρ = Double64[1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]
    # ptranspose 4x4 matrix treated as 2x2 system
    # ptransposing 2nd system (blocks)
    res = ptranspose(ρ, [2, 2], [2])
    @test eltype(res) == Double64
    
    # [1 2; 5 6] (Top left block of input) -> [1 5; 2 6] (Expected top left block of output)
    block1_expected = [1 5; 2 6]
    @test res[1:2, 1:2] == block1_expected
end
end

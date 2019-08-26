@testset "ElementaryArrays" begin
    @testset "kron" begin
        e1 = ElementaryArray{Float64, 2}(1.0, (2, 1), 2, 3)
        e2 = ElementaryArray{Float64, 2}(1.0, (1, 4), 2, 4)
        a = kron(e1, e2)
        b = kron(Array(e1), Array(e2))
        @test a ≈ b
    end
    @testset "matrix addition" begin
        e1 = ElementaryArray{Float64, 2}(1.0, (2, 1), 2, 3)
        e2 = ElementaryArray{Float64, 2}(1.0, (2, 1), 2, 3)
        a = e1 + e2
        b = Array(e1) + Array(e2)
        @test a ≈ b
    end
    @testset "matrix multiplication" begin 
        e1 = ElementaryArray{Float64, 2}(1.0, (2, 1), 2, 3)
        e2 = ElementaryArray{Float64, 2}(1.0, (2, 4), 3, 4)
        a = e1 * e2
        b = Array(e1) * Array(e2)
        @test a ≈ b
    end
end
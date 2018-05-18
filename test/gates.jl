@testset "Random matrices" begin

@testset "QFT" begin
    d = 10
    F = qft(d)
    @test size(F) == (d, d)
    @test norm(F'*F - eye(d)) ≈ 0 atol=1e-13
    @test norm(abs.(F) - fill(1/sqrt(d), d, d)) ≈ 0 atol=1e-13
end

@testset "hadamard" begin
    d = 16
    H = hadamard(d)
    @test size(H) == (d, d)
    @test norm(H'*H - eye(d)) ≈ 0 atol=1e-13
    @test norm(H'-H) ≈ 0 atol=1e-13
    @test_throws ArgumentError hadamard(10)
end

@testset "Pauli matrices" begin
    @test size(sx) == (2, 2)
    @test size(sy) == (2, 2)
    @test size(sz) == (2, 2)

    @test sx*sy - sy*sx == -2im * sz
end

end

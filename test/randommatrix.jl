@testset "Random matrices" begin

@testset "random_dynamical_matrix" begin
    n = 10
    J = zeros(Float64, n^2, n^2)
    random_dynamical_matrix!(J)
    tr = ptrace(J, [n, n], [1])
    @test norm(tr - eye(n)) ≈ 0. atol=1e-13
    @test trace(J) ≈ n atol=1e-13
    @test typeof(J) == Matrix{Float64}
    J = zeros(ComplexF64, n^2, n^2)
    random_dynamical_matrix!(J)
    tr = ptrace(J, [n, n], [1])
    @test norm(tr - eye(n)) ≈ 0. atol=1e-13
    @test trace(J) ≈ n atol=1e-13
    @test typeof(J) == Matrix{ComplexF64}
    J = random_dynamical_matrix(Float64, n)
    tr = ptrace(J, [n, n], [1])
    @test size(J) == (n^2, n^2)
    @test trace(J) ≈ n atol=1e-13
    @test norm(tr - eye(n)) ≈ 0. atol=1e-13
    @test typeof(J) == Matrix{Float64}
    J = random_dynamical_matrix(ComplexF64, n)
    tr = ptrace(J, [n, n], [1])
    @test size(J) == (n^2, n^2)
    @test trace(J) ≈ n atol=1e-13
    @test norm(tr - eye(n)) ≈ 0. atol=1e-13
    @test typeof(J) == Matrix{ComplexF64}
end

@testset "random_unitary" begin
  n = 10
  U = random_unitary(n)
  @test norm(U*U' - eye(n)) ≈ 0 atol=1e-13
end

end

using StatsBase

srand(42)
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

    d1, d2 = 2, 4
    J = random_dynamical_matrix(2, 4)
    @test size(J) == (d1*d2, d1*d2)
    @test trace(J) ≈ d1 atol=1e-13
    tr = ptrace(J, [d2, d1], [1])
    @test norm(tr - eye(d1)) ≈ 0. atol=1e-13
    @test typeof(J) == Matrix{ComplexF64}

    d1, d2 = 2, 4
    J = random_dynamical_matrix(2, 4, 10)
    @test size(J) == (d1*d2, d1*d2)
    @test trace(J) ≈ d1 atol=1e-13
    tr = ptrace(J, [d2, d1], [1])
    @test norm(tr - eye(d1)) ≈ 0. atol=1e-13
    @test typeof(J) == Matrix{ComplexF64}
end

@testset "random_unitary" begin
    n = 10
    U = random_unitary(n)
    @test norm(U*U' - eye(n)) ≈ 0 atol=1e-13

    n = 100
    steps = 100
    r = zeros(steps, n)

    for i=1:steps
        U = random_unitary(n)
        r[i, :] = angle.(eigvals(U))
    end
    r = vec(r)
    h = normalize(fit(Histogram, r, weights(ones(r)), -π:0.1π:π, closed=:left))
    @test all(isapprox.(h.weights, 1/2π, atol=0.01))
end

@testset "random_orthogonal" begin
    n = 10
    O = random_orthogonal(n)
    @test norm(O*O' - eye(n)) ≈ 0 atol=1e-13

    n = 100
    steps = 100
    r = zeros(steps, n)
    for i=1:steps
        U = random_unitary(n)
        r[i, :] = angle.(eigvals(U))
    end
    r = vec(r)
    h = normalize(fit(Histogram, r, weights(ones(r)), -π:0.1π:π, closed=:left))
    @test all(isapprox.(h.weights, 1/2π, atol=0.1))
end

@testset "random GUE/GOE" begin
    n = 10
    @test ishermitian(random_GUE(n))
    @test ishermitian(random_GOE(n))
end

@testset "random_isometry" begin
    n, d = 20, 10
    V = random_isometry(n, d)
    @test norm(V'*V - eye(d)) ≈ 0 atol=1e-13
    @test all([isapprox(v, 0, atol=1e-13) || isapprox(v, 1, atol=1e-13) for v=eigvals(V*V')])
    @test_throws ArgumentError random_isometry(d, n)
end
end

function test_random_dynamical_matrix()
    n = 10
    J = zeros(Float64, n^2, n^2)
    random_dynamical_matrix!(J)
    tr = ptrace(J, [n, n], [1])
    @test_approx_eq_eps norm(tr - eye(n)) 0. 1e-13
    @test_approx_eq_eps trace(J) n 1e-13
    @test typeof(J) == Matrix{Float64}
    J = zeros(Complex128, n^2, n^2)
    random_dynamical_matrix!(J)
    tr = ptrace(J, [n, n], [1])
    @test_approx_eq_eps norm(tr - eye(n)) 0. 1e-13
    @test_approx_eq_eps trace(J) n 1e-13
    @test typeof(J) == Matrix{Complex128}
    J = random_dynamical_matrix(Float64, n)
    tr = ptrace(J, [n, n], [1])
    @test size(J) == (n^2, n^2)
    @test_approx_eq_eps trace(J) n 1e-13
    @test_approx_eq_eps norm(tr - eye(n)) 0. 1e-13
    @test typeof(J) == Matrix{Float64}
    J = random_dynamical_matrix(Complex128, n)
    tr = ptrace(J, [n, n], [1])
    @test size(J) == (n^2, n^2)
    @test_approx_eq_eps trace(J) n 1e-13
    @test_approx_eq_eps norm(tr - eye(n)) 0. 1e-13
    @test typeof(J) == Matrix{Complex128}
end

println("testing random_dynamical_matrix")
test_random_dynamical_matrix()
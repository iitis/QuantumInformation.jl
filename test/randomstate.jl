function test_random_ket()
    ϕ = zeros(Float64, 20)
    ψ = zeros(Complex128, 20)
    random_ket!(ϕ)
    random_ket!(ψ)
    @test length(ϕ) == 20
    @test length(ψ) == 20
    @test typeof(ϕ) == Vector{Float64}
    @test typeof(ψ) == Vector{Complex128}
    @test_approx_eq_eps sum(abs2(ϕ)) 1. 1e-15
    @test_approx_eq_eps sum(abs2(ψ)) 1. 1e-15
    
    ϕ = random_ket(Float64, 100)
    ψ = random_ket(Complex128, 100)
    @test length(ϕ) == 100
    @test length(ψ) == 100
    @test typeof(ϕ) == Vector{Float64}
    @test typeof(ψ) == Vector{Complex128}
    @test_approx_eq_eps sum(abs2(ϕ)) 1. 1e-15
    @test_approx_eq_eps sum(abs2(ψ)) 1. 1e-15
    
    ϕ = random_ket(100)
    @test length(ϕ) == 100
    @test typeof(ϕ) == Vector{Complex128}
    @test_approx_eq_eps sum(abs2(ϕ)) 1. 1e-15
end

function test_random_mixed_state_hs()
    ρ = zeros(Float64, 20, 20)
    σ = zeros(Complex128, 20, 20)
    random_mixed_state_hs!(ρ)
    random_mixed_state_hs!(σ)
    @test size(ρ) == (20, 20)
    @test size(σ) == (20, 20)
    @test typeof(ρ) == Matrix{Float64}
    @test typeof(σ) == Matrix{Complex128}
    @test ishermitian(ρ)
    @test ishermitian(σ)
    @test_approx_eq_eps trace(ρ) 1. 1e-15
    @test_approx_eq_eps trace(σ) 1. 1e-15
    
    ρ = random_mixed_state_hs(Float64, 20)
    σ = random_mixed_state_hs(Complex128, 20)
    @test size(ρ) == (20, 20)
    @test size(σ) == (20, 20)
    @test typeof(ρ) == Matrix{Float64}
    @test typeof(σ) == Matrix{Complex128}
    @test ishermitian(ρ)
    @test ishermitian(σ)
    @test_approx_eq_eps trace(ρ) 1. 1e-15
    @test_approx_eq_eps trace(σ) 1. 1e-15
    
    ρ = random_mixed_state_hs(20)
    @test size(ρ) == (20, 20)
    @test typeof(ρ) == Matrix{Complex128}
    @test ishermitian(ρ)
    @test_approx_eq_eps trace(ρ) 1. 1e-15
end

function test_random_jamiolkowski_state()
    n = 10
    J = zeros(Float64, n^2, n^2)
    random_jamiolkowski_state!(J)
    tr = ptrace(J, [n, n], [1])
    @test_approx_eq_eps norm(tr - eye(n)/n) 0. 1e-13
    @test_approx_eq_eps trace(J) 1. 1e-13
    @test typeof(J) == Matrix{Float64}
    J = zeros(Complex128, n^2, n^2)
    random_jamiolkowski_state!(J)
    tr = ptrace(J, [n, n], [1])
    @test_approx_eq_eps norm(tr - eye(n)/n) 0. 1e-13
    @test_approx_eq_eps trace(J) 1. 1e-13
    @test typeof(J) == Matrix{Complex128}
    J = random_jamiolkowski_state(Float64, n)
    tr = ptrace(J, [n, n], [1])
    @test size(J) == (n^2, n^2)
    @test_approx_eq_eps trace(J) 1. 1e-13
    @test_approx_eq_eps norm(tr - eye(n)/n) 0. 1e-14
    @test typeof(J) == Matrix{Float64}
    J = random_jamiolkowski_state(Complex128, n)
    tr = ptrace(J, [n, n], [1])
    @test size(J) == (n^2, n^2)
    @test_approx_eq_eps trace(J) 1. 1e-13
    @test_approx_eq_eps norm(tr - eye(n)/n) 0. 1e-14
    @test typeof(J) == Matrix{Complex128}
end

println("testing random_ket")
test_random_ket()

println("testing random_mixed_state_hs")
test_random_mixed_state_hs()

println("testing random_jamiolkowski_state")
test_random_jamiolkowski_state()
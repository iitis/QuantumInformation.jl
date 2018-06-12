@testset "Functionals" begin

@testset "trace_norm" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    @test trace_norm(ρ) ≈ 1
    @test trace_norm(sx*ρ) ≈ 1
end

@testset "trace_distance" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]

    @test trace_distance(ρ, ρ) ≈ 0 atol=1e-14
    @test trace_distance(ρ, σ) ≈ 3/(10sqrt(2)) atol=1e-15
end

@testset "hs_norm" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    @test hs_norm(ρ) ≈ sqrt(3)/2
    @test hs_norm(sx*ρ) ≈ sqrt(3)/2
end

@testset "hs_distance" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]

    @test hs_distance(ρ, ρ) ≈ 0 atol=1e-14
    @test hs_distance(ρ, σ) ≈ 3/10 atol=1e-15
end

@testset "fidelity_sqrt" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]
    @test fidelity_sqrt(ρ, σ) ≈ real(trace(sqrtm(sqrtm(ρ) * σ * sqrtm(ρ)))) atol=1e-15
    @test_throws ArgumentError fidelity_sqrt(ones(2, 3), ones(2, 2))
end

@testset "fidelity" begin
    ϕ = ket(0, 2)
    ψ = ket(1, 2)
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]
    @test fidelity(ϕ, ϕ) ≈ 1. atol=1e-15
    @test fidelity(ϕ, ψ) ≈ 0. atol=1e-15
    @test fidelity(ρ, ψ) ≈ 0.75 atol=1e-15
    @test fidelity(ϕ, σ) ≈ 0.4 atol=1e-15
    @test fidelity(ρ, σ) ≈ real(trace(sqrtm(sqrtm(ρ) * σ * sqrtm(ρ))))^2 atol=1e-15
end

@testset "gate_fidelity" begin
    H = [1 1; 1 -1+0im]/sqrt(2)
    @test gate_fidelity(sx, sy) ≈ 0
    @test gate_fidelity(sx, H) ≈ 1/sqrt(2)
    @test gate_fidelity(sy, H) ≈ 0
    @test gate_fidelity(sz, H) ≈ 1/sqrt(2)
    @test gate_fidelity(H, sx) ≈ 1/sqrt(2)
end

@testset "shannon_entropy" begin
    @test shannon_entropy(1/2) ≈ log(2) atol=1e-15
    @test shannon_entropy(1/4) ≈ 0.5623351446188083 atol=1e-15
    @test shannon_entropy(0.5*ones(20)) ≈ 10log(2) atol=1e-15
end

@testset "quantum_entropy" begin
    ϕ = ket(0, 2)
    ρ = [0.25 0.25im; -0.25im 0.75]
    @test quantum_entropy(ϕ) == 0
    @test quantum_entropy(ρ) ≈ 0.25(log(64)-2sqrt(2)*acoth(sqrt(2)))
end

@testset "relative_entropy, js_divergence" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]
    expected = -acoth(5/sqrt(2))/sqrt(2)+acoth(sqrt(2))/sqrt(2)+log(5)-log(46)/2
    expected2 = -1/5*sqrt(2)*atanh(3/(4sqrt(2)))+1/2*log(46/25)
    @test relative_entropy(ρ, σ) ≈ expected
    @test kl_divergence(ρ, σ) ≈ expected
    @test relative_entropy(σ, ρ) ≈ expected2
    @test kl_divergence(σ, ρ) ≈ expected2
    @test js_divergence(ρ, σ) ≈ 0.5expected+0.5expected2
end

@testset "bures_distance, bures_angle" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]
    fsqrt = 1/2 * sqrt(1/5 * (12 + sqrt(46)))

    @test bures_distance(ρ, σ) ≈ sqrt(2 - 2fsqrt)
    @test bures_distance(σ, ρ) ≈ sqrt(2 - 2fsqrt)
    @test bures_angle(ρ, σ) ≈ acos(fsqrt)
    @test bures_angle(σ, ρ) ≈ acos(fsqrt)
end

@testset "superfidelity" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]
    expected = 1/20 * (12 + sqrt(46))
    @test superfidelity(ρ, σ) ≈ expected
    @test superfidelity(ρ, σ) ≈ superfidelity(σ, ρ)
end

@testset "(log)_negativity, ppt" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]
    @test negativity(ρ ⊗ σ, [2, 2], 2) ≈ 0 atol=1e-15
    @test log_negativity(ρ ⊗ σ, [2, 2], 2) ≈ 0 atol=1e-15
    @test ppt(ρ ⊗ σ, [2, 2], 2) > 0

    ρ = proj([1, 0, 0, 1.+0im])/2
    @test negativity(ρ, [2, 2], 2) ≈ 1/2 atol=1e-15
    @test log_negativity(ρ, [2, 2], 2) ≈ log(2) atol=1e-15
    @test ppt(ρ, [2, 2], 2) ≈ -1/2 atol=1e-15

    ρ = 9/10*proj([1, 0, 0, 1.+0im])/2 + 1/10*eye(4)/4

    @test negativity(ρ, [2, 2], 2) ≈ 17/40 atol=1e-15
    @test log_negativity(ρ, [2, 2], 2) ≈ log(37/20) atol=1e-15
    @test ppt(ρ, [2, 2], 2) ≈ -17/40 atol=1e-15
end

@testset "diamond_norm" begin
    srand(42)
    J = random_dynamical_matrix(3, 4)
    @test diamond_norm(J, 3, 4) ≈ 1 atol=1e-6
end

@testset "diamond_distance" begin
    d = 4
    U1 = qft(d)
    U2 = eye(ComplexF64, d)
    J1 = cat([1, 2], [proj(U1[:, i]) for i=1:d]...)
    J2 = cat([1, 2], [proj(U2[:, i]) for i=1:d]...)
    @test diamond_distance(J1, J2, d) ≈ 2 atol=1e-6
end

end

@testset "Utility functions" begin

@testset "renormalize" begin
    v = randn(10)
    renormalize!(v)

    @test norm(v) ≈ 1 atol=1e-13

    A = randn(10, 10)
    renormalize!(A)

    @test trace(A) ≈ 1 atol=1e-13
end

@testset "funcmh" begin
    ρ = [0.25 0.25im; -0.25im 0.75]
    R = zeros(ρ)
    ref = expm(ρ)
    funcmh!(x->exp(x), Hermitian(ρ), R)
    @test R ≈ ref atol=1e-15
    ρ = [0.25 0.25im; -0.25im 0.75]
    R = funcmh(x->exp(x), ρ)
    @test R ≈ ref atol=1e-15

    ρ = [0.25 0.25im; -0.25im 0.75]
    R = zeros(ρ)
    funcmh!(x->exp(x), ρ, R)

    @test R ≈ ref atol=1e-15
end

@testset "random_ball" begin
    b = random_ball(10)

    @test norm(b) < 1.
end

end

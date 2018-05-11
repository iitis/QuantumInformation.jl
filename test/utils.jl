function test_funcmh()
    ρ = [0.25 0.25im; -0.25im 0.75]
    R = zeros(ρ)
    ref = expm(ρ)
    funcmh!(x->exp(x), Hermitian(ρ), R)
    @test R ≈ ref atol=1e-15
    ρ = [0.25 0.25im; -0.25im 0.75]
    R = funcmh(x->exp(x), ρ)
    @test R ≈ ref atol=1e-15
end

println("testing funcmh")
# test_funcmh() TODO: rewrite this function

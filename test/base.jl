function test_ket()
    ϕ = ket(0, 4)
    ψ = ComplexF64[1, 0, 0, 0]
    @test norm(ϕ - ψ) ≈ 0.
    # TODO : Fix these tests: types depend on julia version
    # @test typeof(ket(Float64, 0, 4)) == Vector{Float64}
    # @test typeof(ket(ComplexF64, 0, 4)) == Vector{ComplexF64}
end

function test_bra()
    ϕ = bra(0, 4)
    ψ = ComplexF64[1 0 0 0]
    @test norm(ϕ - ψ) ≈ 0.
    #TODO : Fix these tests: types depend on julia version
    # @test typeof(bra(Float64, 0, 4)) == LinearAlgebra.Adjoint{Float64,Array{Float64,1}}
    # @test typeof(bra(ComplexF64, 0, 4)) == LinearAlgebra.Adjoint{Complex{Float64},Array{Complex{Float64},1}}
end

function test_ketbra()
    ϕψ = ketbra(0, 0, 4)
    αβ = zeros(ComplexF64, 4, 4)
    αβ[1, 1] = 1
    @test norm(ϕψ - αβ) ≈ 0.
    @test typeof(ketbra(Float64, 0, 0, 4)) == Matrix{Float64}
    @test typeof(ketbra(ComplexF64, 0, 0, 4)) == Matrix{ComplexF64}
end

function test_proj()
    ϕ = ket(0, 4)
    ϕϕ = proj(ϕ)
    ψψ = zeros(ComplexF64, 4 ,4)
    ψψ[1, 1] = 1
    @test norm(ϕϕ - ψψ) ≈ 0.
end

function test_base_matrices()
    d = 4
    m = collect(Matrix{ComplexF64}, base_matrices(4))
    for i=1:d, j=1:d
        v = trace(m[i]' * m[j])
        i == j ? @test(v == 1.) : @test(v == 0.)
    end
end

function test_res()
    ρ = [0.25 0.25im; -0.25im 0.75]
    ϕ = res(ρ)
    ψ = [0.25, 0.25im, -0.25im, 0.75]
    @test norm(ϕ - ψ) ≈ 0.
end

function test_unres()
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = unres(res(ρ))
    @test norm(ρ - σ) ≈ 0.
    a = [1 2.1 3; 4 5 6]
    @test norm(unres([1, 2.1, 3, 4 ,5, 6], 2, 3) - a) ≈ 0.
end

function test_kraus_to_superoperator()
    α = 0.25
    K₁ = ComplexF64[0 sqrt(α); 0 0]
    K₂ = ComplexF64[1 0; 0 sqrt(1 - α)]
    kl = Matrix{ComplexF64}[K₁, K₂]
    M = kraus_to_superoperator(kl)
    T = diagm(ComplexF64[1, sqrt(1 - α), sqrt(1 - α), 1 - α])
    T[2, 3] = α
    @test norm(T-M) ≈ 0 atol=1e-15
end

function test_channel_to_superoperator()
    ρ = [0.25 0.25im; -0.25im 0.75]
    T = hcat([ComplexF64[0.25, 0.25im, -0.25im, 0.75] for i=1:4]...) #stack res ρ
    M = channel_to_superoperator(x -> ρ, 2)
    @test norm(T-M) ≈ 0. atol=1e-15
end

function test_apply_kraus()
    α = 0.25
    K₁ = ComplexF64[0 sqrt(α); 0 0]
    K₂ = ComplexF64[1 0; 0 sqrt(1 - α)]
    kl = Matrix{ComplexF64}[K₁, K₂]
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = apply_kraus(kl, ρ)
    ξ = ComplexF64[1/4 + 3/16 sqrt(3/4)*1im/4; -sqrt(3/4)*1im/4 9/16]
    @test trace(σ) ≈ 1. atol=1e-15
    @test ishermitian(σ)
    @test norm(σ - ξ) ≈ 0. atol=1e-15
end

function test_ptrace()
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]
    ξ = ptrace(ρ ⊗ σ, [2, 2], [2,])
    @test norm(ρ - ξ) ≈ 0. atol=1e-15

    ϕ = 1/sqrt(2) * (ket(0, 4) + ket(3, 4))
    ξ = ptrace(proj(ϕ), [2, 2], [2,])
    @test norm(ξ - eye(2)/2) ≈ 0. atol=1e-15
    ξ = ptrace(ϕ, [2, 2], 2)
    @test norm(ξ - eye(2)/2) ≈ 0. atol=1e-15
end

function test_ptranspose()
  ρ =  ComplexF64[1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]
  trans1 = [1 2 9 10; 5 6 13 14; 3 4 11 12; 7 8 15 16]
  trans2 = [1 5 3 7; 2 6 4 8; 9 13 11 15; 10 14 12 16]
  @test norm(ptranspose(ρ, [2, 2], [1]) - trans1) ≈ 0. atol=1e-15
  @test norm(ptranspose(ρ, [2, 2], [2]) - trans2) ≈ 0. atol=1e-15
end

#function test_number2mixedradix()
#    number = 486
#    bases = Int64[8, 42, 2]
#    println(number2mixedradix(number, bases))
#    @test number2mixedradix(number, bases) == Int64[3, 7, 0]
#end

#function test_mixedradix2number()
#    number = Int64[3, 7, 0]
#    bases = Int64[8, 42, 2]
#    @test mixedradix2number(number, bases) == 486
#end

function test_reshuffle()
    X = reshape([1:16;], 4, 4)'
    T = [1 2 5 6; 3 4 7 8; 9 10 13 14; 11 12 15 16]
    @test reshuffle(X) == T
end

function test_trace_distance()
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]

    @test trace_distance(ρ, ρ) ≈ 0 atol=1e-15
    @test trace_distance(ρ, σ) ≈ 0.42426406871192857 atol=1e-15
end

function test_fidelity_sqrt()
    ρ = [0.25 0.25im; -0.25im 0.75]
    σ = [0.4 0.1im; -0.1im 0.6]
    @test fidelity_sqrt(ρ, σ) ≈ real(trace(sqrtm(sqrtm(ρ) * σ * sqrtm(ρ)))) atol=1e-15
end

function test_fidelity()
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

function test_shannon_entropy()
    @test shannon_entropy(1/2) ≈ log(2) atol=1e-15
    @test shannon_entropy(1/4) ≈ 0.562335 atol=1e-15
    @test shannon_entropy(0.5*ones(20)) ≈ 10log(2) atol=1e-15
end

function test_entropy()
    ϕ = ket(0, 2)
    ρ = [0.25 0.25im; -0.25im 0.75]
    @test entropy(ϕ) == 0
    #TODO: add tests for mixed states
end

println("testing ket")
test_ket()

println("testing bra")
test_bra()

println("testing ketbra")
test_ketbra()

println("testing proj")
test_proj()

println("testing base_matrices")
test_base_matrices()

println("testing res")
test_res()

println("testing unres")
test_unres()

println("testing kraus_to_superoperator")
test_kraus_to_superoperator()

println("testing channel_to_superoperator")
test_channel_to_superoperator()

println("testing apply_kraus")
test_apply_kraus()

println("testing ptrace")
test_ptrace()

println("testing ptranspose")
test_ptranspose()

#println("testing number2mixedradix")
#test_number2mixedradix()

#println("testing mixedradix2number")
#test_mixedradix2number()

println("testing reshuffle")
test_reshuffle()

println("testing trace_distance")
test_trace_distance()

println("testing fidelity_sqrt")
test_fidelity_sqrt()

println("testing fidelity")
test_fidelity()

println("testing entropy")
test_entropy()

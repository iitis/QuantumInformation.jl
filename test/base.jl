function test_ket()
    ϕ = ket(0, 4)
    ψ = Complex128[1, 0, 0, 0]
    norm(ϕ - ψ)
end

function test_bra()
    ϕ = bra(0, 4)
    ψ = Complex128[1 0 0 0]
    norm(ϕ - ψ)
end

println("testing ket")
res = test_ket()
println("residuum: $res\n")
@test res < 1e-12

println("testing bra")
res = test_bra()
println("residuum: $res\n")
@test res < 1e-12
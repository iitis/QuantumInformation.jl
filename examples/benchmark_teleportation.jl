using QuantumInformation
using LinearAlgebra
using Statistics

function teleport(d::Int)
    haar = HaarKet(2)
    ψ = reshuffle(𝕀(d))/sqrt(d)
    ϕ = rand(haar)
    had = UnitaryChannel{Matrix{ComplexF64}}(hadamard(d))
    post = [PostSelectionMeasurement(proj(ket(i, d^2)) ⊗ 𝕀(d)) for i in 1:(d ^ 2)]
    ξ = ϕ ⊗ ψ
    ρ =
        ((had ⊗ IdentityChannel(4))∘(cnot ⊗ IdentityChannel(2))∘(IdentityChannel(4) ⊗ Φ(γ)))(
            ξ,
        )
    for j in 1:d
        σ = rots[j](ptrace(post[j](ρ), [2, 2, 2], [1, 2]))
        r[i, k, j] = fidelity(ϕ, σ/tr(σ))
    end
end

rots = [UnitaryChannel(𝕀(2)), UnitaryChannel(sx), UnitaryChannel(sz), UnitaryChannel(sx*sz)]
cnot = UnitaryChannel{Matrix{ComplexF64}}([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0])
r = zeros(steps, length(γs), 4)
mean(r, dims=1)

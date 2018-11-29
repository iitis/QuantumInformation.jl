using QuantumInformation
using LinearAlgebra
using Statistics

steps = 100
haar = HaarKet(2)
ψ = (ket(1, 4) + ket(4, 4))/sqrt(2)
γs = 0.0:0.01:1.0
Φ(γ) = KrausOperators([[1 0; 0 sqrt(1-γ)], [0 sqrt(γ); 0 0]])
post = [PostSelectionMeasurement(proj(ket(i, 4)) ⊗ 𝕀(2)) for i=1:4]
rots = [UnitaryChannel(𝕀(2)), UnitaryChannel(sx), UnitaryChannel(sz), UnitaryChannel(sx*sz)]
had = UnitaryChannel{Matrix{ComplexF64}}(hadamard(2))
cnot = UnitaryChannel{Matrix{ComplexF64}}([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0])
r = zeros(steps, length(γs), 4)
for (k, γ) in enumerate(γs)
    for i = 1:steps
        ϕ = rand(haar)
        ξ = ϕ ⊗ ψ
        ρ = ((had ⊗ IdentityChannel(4))∘(cnot ⊗ IdentityChannel(2))∘(IdentityChannel(4) ⊗ Φ(γ)))(ξ)
        for j=1:4
            σ = rots[j](ptrace(post[j](ρ), [2, 2, 2], [1, 2]))
            r[i, k, j] = fidelity(ϕ, σ/tr(σ))
        end
    end
end
mean(r, dims=1)

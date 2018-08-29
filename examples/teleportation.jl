using QI

steps = 100
haar = HaarKet(2)
ψ = (ket(0, 4) + ket(3, 4))/sqrt(2)
γs = 0.0:0.01:1.0
Φ = KrausOperators([[1 0; 0 sqrt(1-γ)], [0 sqrt(γ); 0 0]])
post = [PostSelectionMeasurement(proj(ket(i, 4)) ⊗ 𝕀(2)) for i=0:3]
rots = [UnitaryChannel(𝕀(2)), UnitaryChannel(sx), UnitaryChannel(sz), UnitaryChannel(sx*sz)]
had = UnitaryChannel{Matrix{Complex128}}(hadamard(2))
cnot = UnitaryChannel{Matrix{Complex128}}([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0])
r = zeros(steps, length(γs), 4);
for (k, γ) in enumerate(γs)
    for i=1:steps
        ϕ = rand(haar)
        ξ = ϕ ⊗ ψ
        ρ = ((had ⊗ IdentityChannel(4))∘(cnot ⊗ IdentityChannel(2))∘(IdentityChannel(4) ⊗ Φ))(ξ)
        for j=1:4
            σ = rots[j](ptrace(post[j](ρ), [2, 2, 2], [1, 2]))
            r[i, k, j] = fidelity(ϕ, σ/trace(σ))
        end
    end
end
mean(r, 1)

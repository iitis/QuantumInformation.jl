```@setup QuantumInformation
using QuantumInformation
```

# Examples
As an example we provide the teleportation protocol in the
presence of noise. Imagine we have an entangled pair of particles in the state
\begin{equation}
|\psi\rangle = \frac{1}{\sqrt{2}} \left( |00\rangle + |11\rangle \right).
\end{equation}

One of the particles stays with Alice and another is sent through a noisy
channel to Bob. As a noise model we chose the amplitude damping channel given by
the Kraus operators
\begin{equation}
K_0 = \begin{matrix}
1 & 0 \\
0 & \sqrt{1-\gamma}
\end{matrix} \quad
K_1 = \begin{matrix}
0 & \sqrt{\gamma} \\
0 & 0
\end{matrix}
.
\end{equation}
The channel has one parameter $\gamma \in [0, 1]$ modeling the strength of the
noise. Assume that Alice possesses a random pure state $|\phi\rangle$ that she
teleports to Bob.

Our examples show the fidelity of the final state at Bob's site averaged over
$100$ random pure initial states. We also check how the parameter $\gamma$
influences this fidelity.

```@repl QuantumInformation
steps = 100
haar = HaarKet(2)
ψ = (ket(0, 4) + ket(3, 4))/sqrt(2)
γs = 0.0:0.01:1.0
Φ = KrausOperators([[1 0; 0 sqrt(1-γ)], [0 sqrt(γ); 0 0]])
post = [PostSelectionMeasurement(proj(ket(i, 4)) ⊗ eye(2)) for i=0:3]
rots = [UnitaryChannel(eye(2)), UnitaryChannel(sx), UnitaryChannel(sz),
UnitaryChannel(sx*sz)]
had = UnitaryChannel{Matrix{ComplexF64}}(hadamard(2))
cnot = UnitaryChannel{Matrix{ComplexF64}}([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0])
r = zeros(steps, length(γs), 4);
for (k, γ) in enumerate(γs)
    for i=1:steps
        ϕ = rand(haar)
        ξ = ϕ ⊗ ψ
        ρ = ((had ⊗ IdentityChannel(4))∘(cnot ⊗
IdentityChannel(2))∘(IdentityChannel(4) ⊗ Φ))(ξ)
        for j=1:4
            σ = rots[j](ptrace(post[j](ρ), [2, 2, 2], [1, 2]))
            r[i, k, j] = fidelity(ϕ, σ/tr(σ))
        end
    end
end
mean(r, 1)
```

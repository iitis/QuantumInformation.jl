```@setup QI
Pkg.add("QI")
importall QI
```

# Functionals

Let us \$\\rho, \\sigma \\in L(\\mathcal{X})\$. Trace norm is defined as \$\\|\\rho\\|_1 = \\mathrm{Tr} \\sqrt{\\rho\\rho^\\dagger}\$ and trace distance is defined based on the trace norm as \$D_1(\\rho,\\sigma)=\\frac{1}{2}\\|\\rho-\\sigma\\|\_1$.

```@repl QI
ψ=(1/sqrt(2)) * (ket(0,2) + ket(1,2))
ϕ=(1/2) * ket(0,2) + (sqrt(3)/2) * ket(1,2)

ρ=proj(ψ)
σ=proj(ϕ)

norm_trace(ρ)
trace_distance(ρ, σ)
```
Now introduce Hilbert–Schmidt norm and distance by \$\\|\\rho\\|_{HS}=\\mathrm{Tr}\\rho^\\dagger \\rho\$ and \$D_{HS}(\\rho,\\sigma)=\\frac{1}{2}\\|\\rho-\\sigma\\|\_{HS}$, respectively.
```@repl QI
norm_hs(ρ)
hs_distance(ρ, σ)
```

Fidelity is a measure of distance of quantum states. It is an example of a
distance measure which is not a metric on the space of quantum states. The
fidelity of two quantum states \$\\rho, \\sigma \in L(\\mathcal{X})\$ is given by
\$F(\\rho,\\sigma)=\\|\\sqrt{\\rho}\\sqrt{\\sigma}\\|_1\$
```@repl QI
fidelity_sqrt(ρ, σ)
fidelity(ρ, σ)
fidelity(ψ, σ)
fidelity(ρ, ϕ)
fidelity(ψ, ϕ)
```

Superfidelity is an upper bound on the fidelity of two quantum states
It is defined by
\$G(\\rho, \\sigma) = \mathrm{Tr}\\rho \\sigma + \\sqrt{1 - \mathrm{Tr}\\rho^2} \\sqrt{1-\mathrm{Tr} \\sigma^2}\$.

```@repl QI
superfidelity(ρ, σ)
```

In order to introduce the diamond norm, we first introduce the notion of the
induced trace norm. Given \$\\Phi \\in T(\\mathcal{X}, \\mathcal{Y})\$ we define its induced trace
norm as \$\\| \\Phi \\|_1 = \\mathrm{max} \\left\\{ \\| \\Phi(X) \\|_1: X \\in L(\\mathcal{X}), \\| X \\|_1 \\leq 1
\\right\\}\$.
The diamond norm of \$\\Phi\$ is defined as
\$
\\| \\Phi \\|_\\diamond = \\| \\Phi \\otimes \\mathbb{I} \\|_1
\$
One important property of the diamond norm is that for Hermiticity-preserving
\$\\Phi \\in T(\\mathcal{X}, \\mathcal{Y})\$ we obtain
\$.
\\| \\Phi \\|_\\diamond = \\max \\left\\{ \\left\\| (\\Phi \\otimes \\mathbb{I})
\\left(|\\psi\\rangle\\langle\\psi| \\right )\\right\\|_1: |\\psi\\rangle \\in \\mathcal{X} \\otimes \\mathcal{Y},
\\langle\\psi|\\psi\\rangle=1 \\right\\}\$.

```@repl QI
K0 = Matrix([1 0; 0 sqrt(1-γ)])
K1 = Matrix([0 sqrt(γ); 0 0])

Φ = KrausOperators([K0,K1])

L0 = Matrix([1 0; 0 sqrt(1-γ)])
L1 = Matrix([0 0; 0 sqrt(γ)])

Ψ = KrausOperators([K0,K1])

norm_diamond(Φ)
diamond_distance(Φ, Ψ)
```

[*Shannon entropy*](https://en.wikipedia.org/wiki/Entropy_(information_theory)) is defined as \$H(\\mathrm{p})=-\\sum_{i=1}^n p_i\\log_2 p_i\$, where
\$\\mathrm{p}=[p_1,...,p_n]\$ is a vector of probabilities.

```@repl QI
p = vec([0.3 0.2 05])

shannon_entropy(p)
shannon_entropy(0.5)
```

For a quantum system described by a state \$\\rho\$, the [*von Neumann entropy*](https://en.wikipedia.org/wiki/Von_Neumann_entropy) is \$S(\\rho)=-\\mathrm{tr} \\rho \\log \\rho\$.
Let \$\\lambda_i\$,  \$0\\geq i < n\$ be eigenvalues of \$\\rho\$, then \$S(\\rho)\$ can be written as \$S(\\rho)=-\\sum_{i=0}^{n-1} \\lambda_i \\log \\lambda_i\$.
```@repl QI
quantum_entropy(0.4 * ρ + 0.6 * σ)
```

One of the measure of distinguishability between two quantum states is a [*qauntum relative entropy*](https://en.wikipedia.org/wiki/Quantum_relative_entropy), called also Kullback–Leibler divergence, defined as
\$S(\\rho\\|\\sigma)=-\\mathrm{Tr}\\rho\\log\\sigma + \\mathrm{Tr}\\rho\\log\\rho\$
```@repl QI
relative_entropy(ρ, σ)
kl_divergence(ρ, σ)
```

Another type of measure of distinguishability between two quantum state is [*quantum Jensen–Shannon divergence*](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence#Quantum_Jensen%E2%80%93Shannon_divergence) given by
\$QJS(\\rho,\\sigma)=S\\left(\\frac{1}{2}\\rho+\\frac{1}{2}\\sigma\\right)-\\left(\\frac{1}{2}S(\\rho)+\\frac{1}{2}S(\\sigma)\\right)\$.
```@repl QI
js_divergence(ρ, σ)
```

[*The Bures distance*](https://en.wikipedia.org/wiki/Bures_metric) defines an infinitesimal distance between quantum states, and it is defined as \$D_B=\\sqrt{2(1-\\sqrt{F(\\rho,\\sigma)})}\$. The value related with Bures distance is Bures angle \$D_A(\\rho,\\sigma)=\\arccos(\\sqrt{F(\\rho,\\sigma)})\$
```@repl QI
bures_distance(ρ, σ)
bures_angle(ρ, σ)
```

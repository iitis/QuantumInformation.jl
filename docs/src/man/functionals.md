```@setup QuantumInformation
using QuantumInformation
using LinearAlgebra
```

# Functionals

### Trace norm and distance
Let $\rho, \sigma \in
\mathrm{L}(\mathcal{X})$. The *trace norm* is defined as $\|\rho\|_1 =
\mathrm{Tr} \sqrt{\rho\rho^\dagger}$ and the trace distance is defined as
$D_1(\rho,\sigma)=\frac{1}{2}\|\rho-\sigma\|_1$.
```@repl QuantumInformation
ψ=(1/sqrt(2)) * (ket(1,2) + ket(2,2))
ϕ=(1/2) * ket(1,2) + (sqrt(3)/2) * ket(2,2)

ρ=proj(ψ)
σ=proj(ϕ)

norm_trace(ρ)
trace_distance(ρ, σ)
```

### Hilbert--Schmidt norm and distance
The [*Hilbert–Schmidt norm*](https://en.wikipedia.org/wiki/Hilbert%E2%80%93Schmidt_operator) norm and distance defined by $\|\rho\|_{HS}=\sqrt{\mathrm{Tr}\rho^\dagger \rho}$ and
$D_{HS}(\rho,\sigma)=\frac{1}{2}\|\rho-\sigma\|_{HS}$, respectively, can be used as follows
```@repl QuantumInformation
norm_hs(ρ)
hs_distance(ρ, σ)
```

### Fidelity and superfidelity
[*Fidelity*](https://en.wikipedia.org/wiki/Fidelity_of_quantum_states) is a measure of distance of quantum states. It is an example of a
distance measure which is not a metric on the space of quantum states. The
fidelity of two quantum states $\rho, \sigma \in \mathrm{L}(\mathcal{X})$ is
given by
$F(\rho,\sigma)=\|\sqrt{\rho}\sqrt{\sigma}\|_1$
```@repl QuantumInformation
fidelity_sqrt(ρ, σ)
fidelity(ρ, σ)
fidelity(ψ, σ)
fidelity(ρ, ϕ)
fidelity(ψ, ϕ)
```
### Superfidelity
[*Superfidelity*](https://www.quantiki.org/wiki/superfidelity) is an upper bound on the fidelity of two quantum states
It is defined by
$G(\rho, \sigma) = \mathrm{Tr}\rho \sigma + \sqrt{1 - \mathrm{Tr}\rho^2}
\sqrt{1-\mathrm{Tr} \sigma^2}$.

```@repl QuantumInformation
superfidelity(ρ, σ)
```

### Diamond norm
In order to introduce the \emph{diamond norm}, we first introduce the notion of
the
induced trace norm. Given $\Phi \in \mathrm{T}(\mathcal{X}, \mathcal{Y})$ we
define its
induced trace
norm as $\| \Phi \|_1 = \mathrm{max} \left\{ \| \Phi(X) \|_1: X \in
L(\mathcal{X}), \| X
\|_1 \leq 1 \right\}$.
The diamond norm of $\Phi$ is defined as
$\| \Phi \|_\diamond = \| \Phi \otimes \mathbb{I}_{\mathrm{L}(\mathcal{Y})} \|_1$
One important property of the diamond norm is that for Hermiticity-preserving
$\Phi \in \mathrm{T}(\mathcal{X}, \mathcal{Y})$ we obtain
$\| \Phi \|_\diamond = \max \left\{ \left\| (\Phi \otimes
\mathbb{I}_{\mathrm{L}(\mathcal{Y})})
	\left(|\psi\rangle\langle\psi| \right )\right\|_1: |\psi\rangle \in
	\mathcal{X} \otimes
	\mathcal{Y},
	\langle\psi|\psi\rangle=1 \right\}$.

```@repl QuantumInformation
γ = 0.3
K0 = Matrix([1 0; 0 sqrt(1-γ)])
K1 = Matrix([0 sqrt(γ); 0 0])

Φ = convert(DynamicalMatrix{Array{Float64,2}}, KrausOperators([K0,K1]))

L0 = Matrix([1 0; 0 sqrt(1-γ)])
L1 = Matrix([0 0; 0 sqrt(γ)])

Ψ = convert(DynamicalMatrix{Array{Float64,2}}, KrausOperators([K0,K1]))

norm_diamond(Φ)
diamond_distance(Φ, Ψ)
```

### Shannon entropy and von Neumann entropy
[*Shannon entropy*](https://en.wikipedia.org/wiki/Entropy_(information_theory)) is defined for a probability vector $p$ as $H(\mathrm{p})=-\sum_{i=1}^n
p_i\log_2 p_i$. We also provide an implementation for the point Shannon entropy.
It is defined as $h(a) = -a \log a - (1-a)\log(1-a)$.

```@repl QuantumInformation
p = vec([0.3 0.2 05])

shannon_entropy(p)
shannon_entropy(0.5)
```

For a quantum system described by a state \$\\rho\$, the [*von Neumann entropy*](https://en.wikipedia.org/wiki/Von_Neumann_entropy) is $S(\rho)=-\mathrm{Tr} \rho \log \rho$. Let $\lambda_i$,  $0\leq i < n$ be the
eigenvalues of $\rho$, then $S(\rho)$ can be written as $S(\rho)=-\sum_{i=1}^{n}
\lambda_i \log \lambda_i$.
```@repl QuantumInformation
ρ = [0.25 0.25im; -0.25im 0.75]
σ = [0.4 0.1im; -0.1im 0.6]

vonneumann_entropy(0.4 * ρ + 0.6 * σ)
```

### Distinguishability between two quantum states
One of the measure of distinguishability between two quantum states is a [*qauntum relative entropy*](https://en.wikipedia.org/wiki/Quantum_relative_entropy), called also Kullback–Leibler divergence, defined as
$S(\rho\|\sigma)=-\mathrm{Tr}\rho\log\sigma + \mathrm{Tr}\rho\log\rho$
```@repl QuantumInformation
relative_entropy(ρ, σ)
kl_divergence(ρ, σ)
```

Another type of measure of distinguishability between two quantum state is [*quantum Jensen–Shannon divergence*](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence#Quantum_Jensen%E2%80%93Shannon_divergence) given by
$QJS(\rho,\sigma)=S\left(\frac{1}{2}\rho+\frac{1}{2}\sigma\right)-
\left(\frac{1}{2}S(\rho)+\frac{1}{2}S(\sigma)\right)$.
```@repl QuantumInformation
js_divergence(ρ, σ)
```

[*The Bures distance*](https://en.wikipedia.org/wiki/Bures_metric) defines an infinitesimal distance between quantum
states, and it is defined as $D_B=\sqrt{2(1-\sqrt{F(\rho,\sigma)})}$. The value
related with Bures distance is the Bures angle
$D_A(\rho,\sigma)=\arccos(\sqrt{F(\rho,\sigma)})$
```@repl QuantumInformation
bures_distance(ρ, σ)
bures_angle(ρ, σ)
```

### Quantum entanglement
One of the entanglement measure is [*negativity*](https://en.wikipedia.org/wiki/Negativity_(quantum_mechanics)) defined as
$\mathrm{N}(\rho)=\frac{\|\rho^{T_A}\|_1-1}{2}$.
```@repl QuantumInformation
negativity(ρ ⊗ σ, [2, 2], 2)
negativity(proj((1/sqrt(2)*(ket(1,2) ⊗ ket(1,2)-ket(2,2) ⊗ ket(2,2)))), [2, 2], 2)

log_negativity(ρ ⊗ σ, [2, 2], 2)
```

[*Positive partial transpose*](https://en.wikipedia.org/wiki/Peres%E2%80%93Horodecki_criterion) (the Peres–Horodecki criterion) is a necessary
condition of separability of the joint state $\rho_{AB}$. According PPT
criterion, if $\rho^{T_B}$ has non negative eigenvalues, then $\rho_{AB}$ is
separable.

```@repl QuantumInformation
ppt(ρ ⊗ σ, [2, 2], 2)
ppt(proj((1/sqrt(2)*(ket(1,2) ⊗ ket(1,2)-ket(2,2) ⊗ ket(2,2)))), [2, 2], 2)
```

Another way to quantification of quantum entanglement is [*Concurrence*](https://en.wikipedia.org/wiki/Concurrence_(quantum_computing)). Concurrence of quantum state $\rho$ is a strong
separability criterion. For two-qubit systems it is defined as
$C(\rho)=\max(0,\lambda_1-\lambda_2-\lambda_3-\lambda_4)$, where $\lambda_i$ are
decreasing eigenvalues of $\sqrt{\sqrt{\rho}\tilde{\rho}\sqrt{\rho}}$ with
$\tilde{\rho}=(\sigma_y\otimes\sigma_y)\rho^*(\sigma_y\otimes\sigma_y)$. If
$C(\rho)=0$, then $\rho$ is separable.
```@repl QuantumInformation
ρ = [0.25 0.1im; -0.1im 0.75]
σ = [0.4 0.1im; -0.1im 0.6]
concurrence(ρ ⊗ σ)
concurrence(proj(max_entangled(4)))
```

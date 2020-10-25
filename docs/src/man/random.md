```@setup QuantumInformation
using QuantumInformation
using LinearAlgebra
```

# Random quantum objects

## Random quantum states

### Pure states
Pure states are elements of the unit sphere in $\mathrm{X}$. Thus it is straightforward
to generate them randomly. We start with a vector of $\dim(\mathrm{X})$ independent
complex numbers sampled from the standard normal distribution. What remains is
to normalize the length of this vector to unity.

This is implemented using the `HaarKet{β}` type. The value $\beta=1$
corresponds to the Haar measure on the unit sphere in $\mathbb{R}^d$, while
$\beta=2$ corresponds to the Haar measure on the unit sphere in $\mathbb{C}^d$. The
usage is as follows
```@repl QuantumInformation
h = HaarKet{2}(3)

ψ = rand(h)

norm(ψ)
```
For convenience we provide the following constructor `HaarKet(d::Int) = HaarKet{2}(d)` as the majority of uses cases require sampling complex states.

### Mixed states
Random mixed states can be generated in one of two equivalent ways. The first
one comes from the partial trace of random pure states. Suppose we have a pure
state $|\psi\rangle \in \mathrm{X} \otimes \mathrm{Y}$. Then we can obtain a random mixed as
\begin{equation}
\rho = \mathrm{Tr}_\mathrm{Y} |\psi\rangle\langle\psi|.
\end{equation}
Note that in the case $\dim(\mathrm{X})=\dim(\mathrm{Y})$ we recover the (flat)
Hilbert-Schmidt distribution on the set of quantum states.

An alternative approach is to start with a Ginibre matrix $G \in \mathrm{L}(\mathrm{X},
\mathrm{Y})$. We obtain a random quantum state $\rho$ as
\begin{equation}
\rho = GG^\dagger/\mathrm{Tr}(GG^\dagger).
\end{equation}
It can be easily verified that this approach is equivalent to the one utilizing
random pure states. First, note that in both cases we start with $\dim(\mathrm{X})
\dim(\mathrm{Y})$ complex random numbers sampled from the standard normal
distribution. Next, we only need to note that taking the partial trace of a
pure state $|\psi\rangle$ is equivalent to calculating $AA^\dagger$ where $A$ is
a matrix obtained from reshaping $|\psi\rangle$.

Sampling random mixed states is implemented using the
`HilbertSchmidtStates{β, K}` type. The meaning of the type parameters
is the same as in the Wishart matrices case. We provide additional constructors
which set the default values of the parameters `HilbertSchmidtStates{β}(d::Int) where β = HilbertSchmidtStates{β, 1}(d)`, `HilbertSchmidtStates(d::Int) = HilbertSchmidtStates{2, 1}(d)`.
The latter one is the most frequent use case. Here is an example
```@repl QuantumInformation
h = HilbertSchmidtStates(3)

ρ = rand(h)

tr(ρ)

eigvals(ρ)
```

### Random quantum channels

Quantum channels are a special subclass of quantum states with constraints
imposed on their *partial* trace as well as trace. Formally, we start with
a Ginibre matrix $G \in \mathrm{L} (\mathrm{X} \otimes \mathrm{Y}, \mathrm{Z})$. We obtain a random
Choi-Jamiołkowski matrix $J_\Phi$ corresponding to a channel $\Phi$ as
$J_\Phi = \left( \mathbb{I}_{\mathrm{X}} \otimes (\mathrm{Tr}_\mathrm{X} GG^\dagger)^{-1/2} \right) GG^\dagger \left( \mathbb{I}_{\mathrm{X}} \otimes (\mathrm{Tr}_\mathrm{X} GG^\dagger)^{-1/2} \right).$

When $\dim(\mathrm{Z})=\dim(\mathrm{X}) \dim(\mathrm{Y})$ this is known to generate a uniform
distribution over the set of quantum
channels.

The implementation uses the type `ChoiJamiolkowskiMatrices{β, K}`. The
parameters $\beta$ and $K$ have the same meaning as in the Wishart matrix case.
Additionally here, the constructor `ChoiJamiolkowskiMatrices{β, K}(idim::Int, odim::Int)  where {β, K}`
takes two parameters--the input and output dimension of the channel. As in the
previous cases we provide some additional constructors for convenience

`function ChoiJamiolkowskiMatrices{β}(idim::Int, odim::Int) where β
    ChoiJamiolkowskiMatrices{β, 1}(idim, odim)
end`,

`function ChoiJamiolkowskiMatrices{β}(d::Int) where β
    ChoiJamiolkowskiMatrices{β}(d, d)
end`,

`function ChoiJamiolkowskiMatrices(idim::Int, odim::Int)
    ChoiJamiolkowskiMatrices{2}(idim, odim)
end`,

`function ChoiJamiolkowskiMatrices(d::Int)
    ChoiJamiolkowskiMatrices(d, d)
end`.

Here is an example of usage
```@repl QuantumInformation
c = ChoiJamiolkowskiMatrices(2, 3)

Φ = rand(c)

ptrace(Φ.matrix, [3, 2],[1])
```
Note that the resulting sample is of type `DynamicalMatrix`.

## [References](@id refs_rand)

[1] B. Collins, I. Nechita, *Random matrix techniques in quantum information theory*, Journal of Mathematical Physics, 2016;57(1):015215.

[2] W. K. Wootters, *Random quantum statesy*, Foundations of Physics, 1990;20(11):1365--1378.

[3] K. Życzkowski, H. J. Sommers, *Induced measures in the space of mixed quantum states*, Journal of Physics A: Mathematical and General, 2001;34(35):7111.

[4] H. J. Sommers, K. Życzkowski, *Statistical properties of random density matrices*, Journal of Physics A: Mathematical and General, 2004;37(35):8457.

[5] Z. Puchała, Ł. Pawela, K. Życzkowski, *Distinguishability of generic quantum states*, Physical Review A, 2016;93(6):062112.

[6] L. Zhang, U. Singh, A. K. Pati, *Average subentropy, coherence and entanglement of random mixed
  quantum states*, Annals of Physics, 2017;377:125--146.

[7] L. Zhang, *Average coherence and its typicality for random mixed quantum states*, Journal of Physics A: Mathematical and Theoretical,
  2017;50(15):155303.

[8] W. Bruzda, V. Cappellini, H. J. Sommers, K. Życzkowski, *Random quantum operations*, Physics Letters A,
  2009;373(3):320--324.

[9] I. Nechita, Z. Puchała, L. Pawela, K. Życzkowski, *Almost all quantum channels are equidistant*, Journal of Mathematical Physics,
  2018;59(5):052201.

[10] L. Zhang, J. Wang, Z. Chen, *Spectral density of mixtures of random density matrices for qubits*, Physics Letters A,
  2018;382(23):1516--1523.

[11] J. Ginibre, *Statistical ensembles of complex, quaternion, and real matrices*, Journal of Mathematical Physics, 1965;6(3):440--449.

[12] M. L. Mehta, *Random matrices. vol. 142.*, Elsevier; 2004.

[13] K. Życzkowski, M. Kuś, *Random unitary matrices*, Journal of Physics A: Mathematical and General, 1994;27(12):4235.

[14] C. Jarlskog, *A recursive parametrization of unitary matrices*, Journal of Mathematical Physics, 2005;46(10):103508.

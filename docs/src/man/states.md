```@setup QuantumInformation
using QuantumInformation
```

# States and channels

In this and the following sections we will denote complex Euclidean spaces
$\mathbb{C}^d$ with $\mathcal{X}$, $\mathcal{Y}$, $\mathcal{Z}$ etc. When needed the dimension of a space $\mathcal{X}$
will be denoted $\mathrm{dim}(\mathcal{X})$. The set of matrices transforming vectors
from $\mathcal{X}$ to $\mathcal{Y}$ will be denoted $\mathrm{L}(\mathcal{X}, \mathcal{Y})$. For simplicity we will
write $\mathrm{L}(\mathcal{X}) \equiv \mathrm{L}(\mathcal{X}, \mathcal{X})$.

## States

By $|\psi\rangle\in\mathcal{X}$ we denote a normed column
vector. Notice that any $|\psi\rangle$ can be expressed as
$|\psi\rangle=\sum_{i=1}^{n} \alpha_i |i\rangle$, where $\sum_{i=1}^{n}
|\alpha_i|^2=1$ and the set $\{|i\rangle\}_{i=1}^{n}$ is the computational
basis.
```@repl QuantumInformation
ket(1,2)
(1/sqrt(2)) * (ket(1,2) + ket(2,2))
```
According to common academic convention, we count the indices of states starting
from one. Following the standard Dirac notation the symbol $\langle\psi|$ denotes
the row vector dual to $|\psi\rangle$. Therefore $|\psi\rangle=\langle\psi|^\dagger$,
where the symbol ${}^\dagger$ denotes the Hermitian conjugation.
```@repl QuantumInformation
bra(2,3)  
```

The inner product of $|\phi\rangle, |\psi\rangle \in \mathcal{X}$ is denoted by
$\langle\psi|\phi\rangle$ and the norm is defined as
$\||\phi\rangle\|=\sqrt{\langle\phi|\phi\rangle}$.
```@repl QuantumInformation
ψ=(1/sqrt(2)) * (ket(1,2) + ket(2,2))
ϕ=(1/2) * ket(1,2) + (sqrt(3)/2) * ket(2,2)

ϕ'*ψ

sqrt(ϕ'*ϕ)
```

The form $|{\psi}\rangle\langle{\phi}|$
denotes outer product of $|{\psi}\rangle$ and $\langle{\phi}|$ from $\mathrm{L}(\mathcal{X})$.
```@repl QuantumInformation
ketbra(2,3,4)
```
Specifically, $|{\psi}\rangle\langle{\psi}|$ is a rank-one projection operator called as *pure state*. Generally, any [*quantum state*](https://en.wikipedia.org/wiki/Qubit) $\rho$ can be expressed as $\rho=\sum_{i=0}^n q_i |i\rangle\langle i|$, where $\sum_{i=0}^n q_i=1$. Notice that $\rho$ is a trace-one positive semi-definite
linear operator *i.e.*: $\rho=\rho^\dagger$, $\rho\geq 0$
and $\mathrm{tr}{\rho}=1$.
```@repl QuantumInformation
proj(ψ)
```
For convenience, the **QuantumInformation.jl** library provides the
implementations of maximally
mixed, maximally entangled and Werner states.
```@repl QuantumInformation
max_entangled(4)
max_mixed(4)
werner_state(4, 0.4)
```

## Non-standard matrix transformations
We will now introduce
reshaping operators, which map matrices to vectors and vice versa. We start with
the mapping $\mathrm{res}:\mathrm{L}(\mathcal{X,Y})\to\mathcal{Y}\otimes\mathcal{X}$, which
transforms the matrix $\rho$ into a vector row by row. More precisely, for
dyadic operators $|\psi\rangle\langle\phi|$, where $|\psi\rangle \in
\mathcal{Y}$, $|\phi\rangle \in \mathcal{X}$ the operation $\mathrm{res}$ is
defined as
$\mathrm{res}(|\psi\rangle\langle\phi|)=|\psi\rangle|\overline{\phi}\rangle$ and
can be uniquely extend to the whole space $\mathrm{L}(\mathcal{X,Y})$ by
linearity.
```@repl QuantumInformation
res(ketbra(1,2,2))
```
The inverse operation to $\mathrm{res}$ is
$\mathrm{unres}:\mathcal{Y}\otimes\mathcal{X}\to \mathrm{L}(\mathcal{X,Y}) $,
which transforms the vector into a matrix. It is defined as the unique linear
mapping satisfying $\rho=\mathrm{unres}(\mathrm{res}(\rho))$.
```@repl QuantumInformation
unres(res(ketbra(1,2,2)))
```
Let us recall that trace is a mapping
$\mathrm{Tr}:\mathrm{L}(\mathcal{X})\to \mathbb{C},$ given by
$\mathrm{Tr}:\rho\mapsto\sum_{i=1}^{\mathrm{dim}(\mathcal{X})}\langle
e_i|\rho|e_i\rangle$, where $\{|e_i\rangle \}$ is an orthonormal basis of
$\mathcal{X}$. According to this, *partial trace* is a mapping
$\mathrm{Tr}_{\mathcal{X}}: \mathrm{L}(\mathcal{X}\otimes\mathcal{Y}) \to
\mathrm{L}(\mathcal{Y})$ such that $\mathrm{Tr}_{\mathcal{X}}: \rho_A\otimes
\rho_B \mapsto \rho_B \mathrm{Tr}(\rho_A)$, where $\rho_A\in
\mathrm{L}(\mathcal{X})$, $\rho_B\in \mathrm{L}(\mathcal{Y})$. As this is a
linear map, it may be uniquely extended to the case of operators which are not
in a tensor product form.
```@repl QuantumInformation
ρ = [0.25 0.25im; -0.25im 0.75]
σ = [0.4 0.1im; -0.1im 0.6]

ptrace(ρ ⊗ σ, [2, 2], [2])
```
Matrix transposition is a mapping ${}^T:\mathrm{L}(\mathcal{X,Y}) \to
\mathrm{L}(\mathcal{Y,X})$ such that $\left(\rho^T \right)_{ij} = \rho_{ji}$,
where $\rho_{ij}$ is a $i$-th row, $j$-th column element of matrix $\rho$.
Following this, we may introduce \emph{partial transposition}
${}^{\Gamma_B}:
\mathrm{L}(\mathcal{X}_A \otimes \mathcal{X}_B, \mathcal{Y}_A \otimes \mathcal{Y}_B) \to
\mathrm{L}(\mathcal{X}_A \otimes \mathcal{Y}_B, \mathcal{Y}_A \otimes \mathcal{X}_B)$,
which for a product
state $\rho_A\otimes\rho_B$ is given by ${}^{\Gamma_B}:
\rho_A\otimes\rho_B\mapsto\rho_A\otimes\rho_B^T$. The definition of partial
transposition can be uniquely extended for all operators from linearity.
```@repl QuantumInformation
ptranspose(ρ ⊗ σ, [2, 2], [1])
```
For given multiindexed matrix $\rho_{(m,\mu),(n,\nu)}=\langle
m \mu|\rho|n \nu\rangle$, the reshuffle operation is defined as
$\rho^R_{(m,\mu),(n,\nu)}=\rho_{(m,n),(\mu,\nu)}$.
```@repl QuantumInformation
reshuffle(ρ ⊗ σ)
```

## Channels

Physical transformations of quantum states into quantum
states are called quantum channels *i.e.* linear Completely Positive
Trace Preserving (CP-TP) transformations. Probabilistic transformations of
quantum states are called quantum operations and mathematically they are defined
as linear Completely Positive Trace Non-increasing (CP-TNI) maps. For the sake
of simplicity we will refer to both CP-TP and CP-TNI maps as quantum channels
when it will not cause confusion.

There exists various representations of quantum channels such as:
* Kraus operators,
* natural representation, also called superoperator representation,
* Stinespring representation,
* Choi-Jamiołkowski matrices, sometimes called dynamical matrices.

The product of superoperators $\Phi_1\in
\mathrm{T}(\mathcal{X}_1,\mathcal{Y}_1)$, $\Phi_2\in
\mathrm{T}(\mathcal{X}_2,\mathcal{Y}_2)$ is a mapping $\Phi_1\otimes\Phi_2\in
T(\mathcal{X}_1\otimes\mathcal{X}_2,\mathcal{Y}_1\otimes\mathcal{Y}_2)$ that
satisfies
$(\Phi_1\otimes\Phi_2)(\rho_1\otimes\rho_2)=\Phi_1(\rho_1)\otimes\Phi_2(\rho_2)$.
For the operators that are not in a tensor product form this notion can be
uniquely extended from linearity.

According to Kraus' theorem, any completely positive trace-preserving (CP-TP)
map $\Phi$ can always be written as $\Phi(\rho)=\sum_{i=1}^r K_i \rho
K_i^\dagger$ for some set of operators $\{K_i\}_{i=1}^r$ satisfying
$\sum_{i=1}^r K_i^\dagger K_i = \mathbb{I}_\mathcal{X}$, where $r$ is the rank
of superoperator $\Phi$.

Another way to represent the quantum channel is based on Choi-Jamiołkowski
isomorphism. Consider mapping $J:\mathrm{T}(\mathcal{X,Y})\to
\mathrm{L}(\mathcal{Y}\otimes\mathcal{X})$ such that
$J(\Phi)=(\Phi\otimes\mathbb{I}_{\mathrm{L}(\mathcal{X})})
(\mathrm{res}(\mathbb{I}_{\mathcal{X}})
\mathrm{res}(\mathbb{I}_{\mathcal{X}})^\dagger)$.
Equivalently
$J(\Phi)=\sum_{i,j=1}^{\mathrm{dim(\mathcal{X})}}\Phi(|i\rangle\langle
j|)\otimes|i\rangle\langle j|$. The action of a superoperator in the Choi
representation is given by
$\Phi(\rho)=\mathrm{Tr}_\mathcal{X}(J(\Phi)(\mathbb{I}_\mathcal{Y}\otimes\rho^T))$.

The natural representation of a quantum channel $\mathrm{T}(\mathcal{X}, \mathcal{Y})$ is a
mapping $\mathrm{res}(\rho) \mapsto \mathrm{res}(\Phi(\rho))$. It is represented
by a matrix $K(\Phi) \in \mathrm{L}(\mathcal{X} \otimes \mathcal{X}, \mathcal{Y} \otimes \mathcal{Y})$ for which the following holds
\begin{equation}
K(\Phi) \mathrm{res}(\rho) = \mathrm{res}(\Phi(\rho)),
\end{equation}
for all $\rho \in \mathrm{L}(\mathcal{X})$.

Let $\mathcal{X}, \mathcal{Y}$ and $\mathcal{Z}$ be a complex Euclidean spaces. The action of the
Stinespring representation of a quantum channel $\Phi\in \mathrm{T}(\mathcal{X},\mathcal{Y})$ on
a state $\rho\in \mathrm{L}(\mathcal{X})$ is given by
\begin{equation}
\Phi(\rho)=\mathrm{Tr}_\mathcal{Z}(A\rho A^\dagger),
\end{equation}
where $A\in\mathrm{L}(\mathcal{X},\mathcal{Y}\otimes\mathcal{Z})$.

We now briefly describe the relationships among channel representations
[[1]](@ref refs_sc). Let $\Phi\in \mathrm{T}(\mathcal{X}, \mathcal{Y})$
be a quantum channel which can be written in the Kraus representation as
$\Phi(\rho)=\sum_{i=1}^r K_i \rho K_i^\dagger$,
where $\{K_i\}_{i=1}^r$ are Kraus operators satisfying $\sum_{i=1}^r K_i^\dagger
K_i = \mathbb{I}_\mathcal{X}$. According to this assumption, $\Phi$ can be represented
in
* Choi representation as $J(\Phi)=\sum_{i=1}^r \mathrm{res}(K_i)\mathrm{res}(K_i^\dagger)$,
* natural representation as $K(\Phi)=\sum_{i=1}^r K_i\otimes K_i^{*}$,
* Stinespring representation as $\Phi(\rho)=\mathrm{Tr}_\mathcal{Z}(A\rho A^\dagger)$,
where $A=\sum_{i=1}^r K_i\otimes |e_i\rangle$ and $\mathcal{Z}=\mathbb{C}^r$.

In **QuantumInformation.jl** states and channels are always represented in the
computational basis therefore channels are stored in the memory as either
vectors of matrices in case of Kraus operators or matrices in other cases. In
**QuantumInformation.jl** quantum channels are represented by a set of types deriving from
an abstract type `AbstractQuantumOperation{T}` where type parameter
`T` should inherit from `AbstractMatrix{<:Number}`. Every
type inheriting from `AbstractQuantumOperation{T}` should contain
fields `idim` and `odim` representing the dimension of input
and output space of the quantum channel.

Two special types of channels are implemented: `UnitaryChannel` and
`IdentityChannel` that can transform ket vectors into ket vectors.


### Constructors
Channel objects can be constructed from matrices that represent them, as shown
in the following listing
```@repl QuantumInformation
γ=0.4
K0 = Matrix([1 0; 0 sqrt(1-γ)])
K1 = Matrix([0 sqrt(γ); 0 0])

Φ = KrausOperators([K0,K1])

iscptp(Φ)
```
There are no checks whether a matrix represents a valid CP-TP or CP-TNI map,
because this kind of verification is costly and requires potentially expensive
numerical computation. Function such as `iscptp()`, and
`iscptni()` are provided to test properties of supposed quantum
channel or quantum operation.

### Conversion
Conversions between all quantum channel types,
*i.e.* these that derive from `AbstractQuantumOperation{T}`
are implemented. The users are not limited by any single channel representation
and can transform between representations they find the most efficient or
suitable for their purpose.

```@repl QuantumInformation
Ψ1 = convert(SuperOperator{Matrix{ComplexF64}}, Φ)
Ψ2 = convert(DynamicalMatrix{Matrix{Float64}}, Φ)
Ψ3 = convert(Stinespring{Matrix{Float64}}, Φ)
```

### Application
Channels can act on pure and mixed states
represented by vectors and matrices respectively. Channels are callable and therefore mimic
application of a function on a quantum state.

```@repl QuantumInformation
ρ1=ψ * ψ'
Φ(ρ1)
Ψ1(ρ1)
Φ(ψ)
```

### Composition
Channels can be composed in parallel or in
sequence. Composition in parallel is done using `kron()` function or
the overloaded $\otimes$ operator. Composition in sequence can be done in two
ways either by using Julia built-in function composition operator
$(f\circ g)(\cdot)=f(g)(\cdot)$ or by using multiplication of objects inheriting
from `AbstractQuantumOperation{T}` abstract type.

```@repl QuantumInformation
ρ2=ϕ * ϕ'

(Φ⊗Φ)(ρ1⊗ρ2)
(Ψ1∘Ψ2)(ρ1)
```
## [References](@id refs_sc)

[1] J. Watrous, *The Theory of Quantum Information*, Cambridge University Press (2018).

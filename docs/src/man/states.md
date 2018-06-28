```@setup QI
Pkg.add("QI")
importall QI
```

# States and channels

## Generalized transpositions

## States
Let \$\\mathcal{X}\$ be a complex Euclidean space with a standard inner product. By
\$|\\psi\\rangle\in\\mathcal{X}\$ we denote a normed column vector. Notice that
any \$|\\psi\\rangle\$ can be expressed as \$|\\psi\\rangle=\\sum_{i=0}^n \\alpha_i |i\\rangle\$, where \$\\sum_{i=0}^n |\\alpha_i|^2=1\$ and set \$\\{|i\\rangle\\}_{i=0}^n\$ is a computational basis.
```@repl QI
ket(0,2)  
(1/sqrt(2)) * (ket(0,2) + ket(1,2))
```
By \$\\langle\\psi|\$ the
row vector dual to \$|\\psi\\rangle\$ is denoted. It means that \$|\\psi\\rangle=\\langle\\psi|^\\dagger\$, where the symbol \${}^\\dagger\$ denotes the
Hermitian conjugation.
```@repl QI
bra(1,3)  
```

Through the \$\\langle\\psi|\\phi\\rangle\$ is denoted a inner product and the
norm defined as \$\\|\\phi\\rangle\\|=\\sqrt{\\langle\\psi|\\psi\\rangle}\$.
```@repl QI
ψ=(1/sqrt(2)) * (ket(0,2) + ket(1,2))
ϕ=(1/2) * ket(0,2) + (sqrt(3)/2) * ket(1,2)

ϕ'*ψ

sqrt(ϕ'*ϕ)
```

The form \$|{\\psi}\\rangle\\langle{\\phi}|\$
denotes outer product of \$|{\\psi}\\rangle\$ and \$\\langle{\\phi}|\$ from \$\\mathrm{L}(\\mathcal{X})\$, where  \$\\mathrm{L}(\\mathcal{X})\$ is the set of linear operators
from \$\\mathcal{X}\$ to \$\\mathcal{X}\$. In the case of \$\\mathrm{L}(\\mathcal{X,Y})\$ we have as set of linear operators
from \$\\mathcal{X}\$ to \$\\mathcal{Y}\$ and we assume that \$\\mathrm{L}(\\mathcal{X}) :=\\mathrm{L}(\\mathcal{X,X})\$.
```@repl QI
ketbra(1,3,4)
```
Especially, \$|{\\psi}\\rangle\\langle{\\psi}|\$ is a rank-one projection operator called as *pure states*. Generally, any [*quantum state*](https://en.wikipedia.org/wiki/Qubit) \$\\rho\$ can be expressed as \$\\rho=\\sum_{i=0}^n q_i |i\\rangle\\langle i|\$, where \$\\sum_{i=0}^n q_i=1\$. Notice that \$\rho\$ is a trace-one positive semi-definite
linear operator *i.e.*: \$\\rho=\\rho^\\dagger\$, \$\\rho\\geq 0\$
and \$\\mathrm{tr}{\\rho}=1\$.
```@repl QI
proj(ψ)
```
One of the useful operation is matrix vectorization \$\\mathrm{vec}\$, which is a transformation
\$\\mathrm{vec}:\\mathrm{L}(\\mathcal{X,Y})\\mapsto\\mathcal{Y}\\otimes\\mathcal{X}\$ such that
\$\\mathrm{vec}(\\rho)=[a_{1,1},...,a_{m,1}.a_{1,2},...,a_{m,2},...,a_{1,n},....,a_{m,n}]^T\$,
where \$\\rho=[\\rho_{i,j}]\$.
```@repl QI
vec(ketbra(0,1,2))
```
In the opposite to \$\\mathrm{vec}\$ transformation, we have reshaping map \$\\mathrm{res}:\\mathrm{L}(\\mathcal{X,Y})\\mapsto\\mathcal{Y}\\otimes\\mathcal{X}\$, which transform matrix
\$\\rho\$ into a vector row by row. More precisely, for dyadic operators \$|\\psi\\rangle\\langle\\phi|\$, where \$|\\psi\\rangle \\in \\mathcal{X}\$, \$|\\phi\\rangle
\\in \\mathcal{Y}\$ operation \$\\mathrm{vec}\$ is defined as \$\\mathrm{vec}(|\\psi\\rangle\\langle\\phi|)=|\\psi\\rangle|\\overline{\\phi}\\rangle\$ and can be uniquely extend the definition to the whole space \$\\mathrm{L}(\\mathcal{X,Y})\$ by
linearity. It is easy to check that \$\\mathrm{res}(\\rho)=\\mathrm{vec}(\\rho^T)\$.
```@repl QI
res(ketbra(0,1,2))
```

Inverse operation to \$\\mathrm{res}\$ is a \$\\mathrm{unres}:\\mathcal{Y}\\otimes\\mathcal{X}\\mapsto \\mathrm{L}(\\mathcal{X,Y}) \$ map, which transforms
the vector into a matrix. It means that \$\\rho=\\mathrm{unres}(\\mathrm{res}(\\rho))\$.
```@repl QI
unres(res(ketbra(0,1,2)))
```
Let us recall that trace is a mapping \$\\mathrm{Tr}:\\mathrm{L}(\\mathcal{X})\\ \\mapsto \\mathbb{C},\$ given by \$\\mathrm{Tr}:\\rho\\mapsto\\sum_{i=1}^{\\mathrm{dim}(\\mathcal{X})}\\langle e_i|\\rho|e_i\\rangle\$, where \$\\{|e_i\\rangle \\}\$ is an orthonormal basis of \$\\mathcal{X}\$. According this, *partial trace* is a mapping \$\\mathrm{Tr}_{\\mathcal{X}}: \\mathrm{L}(\\mathcal{X}\\otimes\\mathcal{Y}) \\mapsto \\mathrm{L}(\\mathcal{Y})\$ such that \$\\mathrm{Tr}_{\\mathcal{X}}: \\rho_A\\otimes \\rho_B \\mapsto \\rho_B \\mathrm{Tr}(\\rho_A)\$, where \$\\rho_A\\in \\mathrm{L}(\\mathcal{X})\$, \$\\rho_B\\in \\mathrm{L}(\\mathcal{Y})\$.
As this is a linear map, it may be uniquely extended to the case of operators
which are not in a tensor product form.
```@repl QI
ρ = [0.25 0.25im; -0.25im 0.75]
σ = [0.4 0.1im; -0.1im 0.6]

ptrace(ρ ⊗ σ, [2, 2], [2,])
```

Matrix transposition is a mapping \${}^T:\\mathrm{L}(\\mathcal{X,Y}) \\mapsto \\mathrm{L}(\\mathcal{Y,X})\$ such that \$\\rho_{ij}^T = \\rho_{ji}\$, where \$\\rho_{ij}\$ is a \$i\$-th row, \$j\$-th column element of matrix \$\\rho\$. Following this, we may introduce *partial transposition* \${}^{T_B}:\\mathrm{L}(\\mathcal{X_A,Y_A}\\otimes\\mathcal{X_B,Y_B}) \\mapsto \\mathrm{L}(\\mathcal{Y_A,X_A}\\otimes\\mathcal{X_B,Y_B})\$,
which for product state \$\\rho_A\\otimes\\rho_B\$ is given by \${}^{T_A}: \\rho_A\\otimes\\rho_B\\mapsto\\rho_A^T\\otimes\\rho_B\$.
The definition of partial trace can be extended for all operators by trace linearity.
```@repl QI
ptranspose(ρ ⊗ σ, [2, 2], [1])
```

For given multiindexed matrix \$\\rho_{(m,\\mu),(n,\\nu)}=\\langle m|\\langle\\mu|\\rho|n\\rangle\\nu\\rangle\$,
reshuffle operation is defined as \$\\rho^R_{(m,\\mu),(n,\\nu)}=\\rho_{(m,n),(\\mu,\\nu)}\$.
```@repl QI
reshuffle(ρ ⊗ σ)
```

## Channels
In the most general case evolution of the quantum system can be described
using the notion of a [*quantum channel*](https://en.wikipedia.org/wiki/Quantum_channel).
First, introduce a *superoperator* as a linear mapping acting on linear operators \$\\mathrm{L}(\\mathcal{X})\$
and transforming them into operators acting on \$\\mathrm{L}(\\mathcal{Y})\$. The set of all such mapping wiil be denoted by \$\\mathrm{T}(\\mathcal{X},\\mathcal{Y})\$ and \$\\mathrm{T}(\\mathcal{X}):=\\mathrm{T}(\\mathcal{X},\\mathcal{X})\$. In mathematical terms,
a quantum channel is a superoperator \$\\Phi:\\mathrm{L}(\\mathcal{X})\\mapsto \\mathrm{L}(\\mathcal{Y})\$
that is *trace-preserving* (\$\\forall \\rho\\in \\mathrm{L}(\\mathcal{X})\\quad \\mathrm{Tr}(\\Phi(\\rho))=\\mathrm{Tr}(\\rho)\$)
and *completely positive* (\$\\forall \\mathcal{Z} \\forall \\rho \\in \\mathrm{L}(\\mathcal{X\\otimes Z}), \\rho\\geq 0, \\quad \\Phi\\otimes\\mathbb{I}_{\\mathrm{L}(\\mathcal{X})}(\\rho)\\geq 0\$).
The the product of the given super-operators \$\\Phi_1\\in \\mathrm{T}(\\mathcal{X_1},\\mathcal{Y_1})\$, \$\\Phi_2\\in \\mathrm{T}(\\mathcal{X_2},\\mathcal{Y_2})\$ is a mapping \$\\Phi_1\\otimes\\Phi_2\\in T(\\mathcal{X}_1\\otimes\\mathcal{X}_2,\\mathcal{Y}_1\\otimes\\mathcal{Y}_2)\$ that satisfies \$(\\Phi_1\\otimes\\Phi_2)(\\rho_1\\otimes\\rho_2)=\\Phi_1(\\rho_1)\\otimes\\Phi_2(\\rho_2)\$.

According to Kraus' theorem, any completely positive trace-preserving (CPTP) map \$\\Phi\$ can always be written as
\$\\Phi(\\rho)=\\sum_{i=0}^rK_i \\rho K_i^\\dagger\$ for some set of operators \$\\{K_i\\}_{i=0}^r\$ satisfying \$\\sum_i K_i^\\dagger K_i = \\mathbb{I}\$, where \$r\$ is the rank of super-operator \$\\Phi\$. Another way to represent the quantum channel is based on Choi-Jamiołkowski isomorphism. Consider mapping \$J:\\mathrm{T}(\\mathcal{X,Y})\\mapsto \\mathrm{L}(\\mathcal{Y}\\otimes\\mathcal{X})\$ such that \$J(\\Phi)=(\\Phi\\otimes\\mathbb{I}_{\\mathrm{L}(\\mathcal{X})})(\\mathrm{res}(\\mathbb{I}_{\\mathcal{X}})\\mathrm{res}(\\mathbb{I}_{\\mathcal{X}})^\\dagger)\$. Equivalently \$J(\\Phi)=\\sum_{i,j=0}^{\\mathrm{dim(\\mathcal{X})-1}}\\Phi(|i\\rangle\\langle j|)\\otimes|i\\rangle\\langle j|\$. The action of a superoperator in the Choi representation,also called dynamical matrix of \$\\Phi\$, is given by \$\\Phi(\\rho)=\\mathrm{Tr}_\\mathcal{X}(J(\\Phi)(\\mathbb{I}_\\mathcal{Y}\\otimes\\rho^T))\$. Last representation of quantum channel implemented in `QI.jl` is  Stinespring representation. Supoose that \$A\\in \\mathrm{L}(\\mathcal{X},\\mathcal{Y}\otimes\\mathcal{Z})\$, then \$\\Phi(\\rho)=\\mathrm{Tr}_\\mathcal{Z}(A\\rho A^\\dagger)\$.

### Constructors
Channel objects can be constructed from matrices that represents them. As shown in the following listing
```@repl QI
γ=0.4
K0 = Matrix([1 0; 0 sqrt(1-γ)])
K1 = Matrix([0 sqrt(γ); 0 0])

Φ = KrausOperators([K0,K1])

iscptp(Φ)
```
### Conversion
Conversions between all quantum channel types are implemented. The user is not
limited by any single channel representation and can transform between
representations he finds the most efficient or suitable for his purpose.

```@repl QI
Ψ1 = convert(SuperOperator{Matrix{ComplexF64}}, Φ)

# or

SuperOperator{Matrix{ComplexF64}}(Φ)

Ψ2 = convert(DynamicalMatrix{Matrix{Float64}}, Φ)

#or

DynamicalMatrix{Matrix{Float64}}(Φ)

Ψ3 = convert(Stinespring{Matrix{Float64}}, Φ)

#or

Stinespring{Matrix{Float64}}(Φ)
```

### Application
Channels can act on pure and mixed states as represented by vectors and matrices. Channels
are callable and therefore mimic application of a~function on a~quantum state.

```@repl QI
ρ1=ψ * ψ'

Φ(ρ1)
Ψ1(ρ1)
Φ(ψ)
```

### Composition
Channels can be composed in parallel or in sequence. Composition in parallel is done using
`kron()` function or overloaded \$\\otimes\$ operator. Composition in sequence can
be done in two ways either by using `Julia` build in function composition operator \$(f\\circ g)(\\cdot)=f(g)(\\cdot)\$ or by using multiplication of objects inheriting from `AbstractQuantumOperation{T}` abstract type.

```@repl QI
ρ2=ϕ * ϕ'

(Φ⊗Φ)(ρ1⊗ρ2)
(Ψ1∘Ψ2)(ρ1)
```

## Measurement
Measurement is modeled in two ways:
* as Positive Operator Valued Measure (POVM) based,
* measurement with post-selection.
In both cases a~measurement is treated as a special case of quantum channel
(operation) respectively as defined below.

### Positive Operator Valued Measure measurement
A POVM measurement is defined as follows, let a
mapping from a finite alphabet of measurement outcomes to the set of linear
positive operators \$\\mu:\\Gamma\\rightarrow\\mathcal{P}(\\mathcal{X})\$
 be given, if \$\\sum_{\\xi\\in\\Gamma} {\\mu(\\xi)=\\mathbb{I}_{\\mathcal{X}}}\$ then \$\\mu\$ is a POVM measurement. POVM
measurement models the situation where a quantum object is destroyed during the
measurement process and quantum state after the measurement does not exists.

We model POVM measurement as a channel
\$\\theta:\\mathcal{T}(\\mathcal{X})\\rightarrow \\mathcal{T}(\\mathcal{Y})\$, where
\$\\mathcal{Y}=\\mathrm{span}\\{|\\xi\\rangle\\}_{\\xi\\in\\Gamma}\$ such that \$\\theta(\\rho)
= \\sum_{\\xi\\in\\Gamma} \\mathrm{tr}(\\rho\\, \\mu(\\xi))\\|\\xi\\rangle\\langle \\xi|\$. This channel transforms
the measured quantum state into a~classical state (diagonal matrix) containing
probabilities of reassuring given outcomes. Note that in `QI`
\$\\Gamma=\\{1,2,\\ldots,|\\Gamma|\\}\$ and POVM measurements are represented by the
type `POVMMeasurement{T} <: AbstractQuantumOperation{T}` where `T<:AbstractMatrix{<:Number}`.
Predicate function `ispovm()` verifies whether a~list of matrices is a proper POVM.
```@repl QI
ρ=proj(1./sqrt(2)*(ket(0,3)+ket(2,3)))

E0 = proj(ket(0,3))
E1 = proj(ket(1,3))+proj(ket(2,3))

M = POVMMeasurement([E0,E1])

ispovm(M)
M(ρ)
```
When a quantum system after being measured is not destroyed one can be
interested by its state after the measurement. This state depends on the
measurement outcome. In this case the measurement process is defined in the following way.

Let $\mu:\Gamma\rightarrow\mathcal{L}(\mathcal{X}, \mathcal{Y})$
be a mapping from a finite set of measurement outcomes to set of linear operators called effects, then
if $\sum_{\xi\in\Gamma} {\mu(\xi)^\dagger \mu(\xi)=\mathbb{I}_{\mathcal{X}}}$
then $\mu$ is a quantum measurement. Given outcome $\xi$ was obtained, the state before the measurement $\rho$
is transformed into sub-normalized quantum state $\rho_\xi=\mu(\xi)\rho\mu(\xi)^\dagger$. The outcome $\xi$ will be obtained
with probability $\mathrm{Tr}(\rho_\xi)$.
```@repl QI
PM = PostSelectionMeasurement(E1)
iseffect(PM)
PM(ρ)
```
In `QI` this kind of measurement is modeled as CP-TNI map with single Kraus operator $\mu(\xi)$ and represented as
`PostSelectionMeasurement{T} <: AbstractQuantumOperation{T}` where `T<:AbstractMatrix{<:Number}`. Measurement types can be composed and converted to Kraus operators,
superoperators, Stinespring representation operators, and dynamical matrices.
```@repl QI
α = 0.3
K0 = ComplexF64[0 0 sqrt(α); 0 1 0; 0 0 0]
K1 = ComplexF64[1 0 0; 0 0 0; 0 0 sqrt(1 - α)]
Φ = KrausOperators([K0,K1])

ρ=proj(1./sqrt(2)*(ket(0,3)+ket(2,3)))

(PM∘Φ)(ρ)
```

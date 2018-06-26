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
The form \$|{\\psi}\\rangle\\langle{\\phi}|\$
denotes outer product of \$|{\\psi}\\rangle\$ and \$\\langle{\\phi}|\$ from \$\\mathcal{L}(\\mathcal{X})\$, where  \$\\mathcal{L}(\\mathcal{X})\$ is the set of linear operators
from \$\\mathcal{X}\$ to \$\\mathcal{X}\$. In the case of \$\\mathcal{L}(\\mathcal{X,Y})\$ we have as set of linear operators
from \$\\mathcal{X}\$ to \$\\mathcal{Y}\$ and we assume that \$\\mathcal{L}(\\mathcal{X}) :=\\mathcal{L}(\\mathcal{X,X})\$.
```@repl QI
ketbra(1,3,4)
```
Especially, \$|{\\psi}\\rangle\\langle{\\psi}|\$ is a rank-one projection operator called as *pure states*. Generally, any [*quantum state*](https://en.wikipedia.org/wiki/Qubit) \$\\rho\$ can be expressed as \$\\rho=\\sum_{i=0}^n q_i |i\\rangle\\langle i|\$, where \$\\sum_{i=0}^n q_i=1\$. Notice that \$\rho\$ is a trace-one positive semi-definite
linear operator *i.e.*: \$\\rho=\\rho^\\dagger\$, \$\\rho\\geq 0\$
and \$\\mathrm{tr}{\\rho}=1\$.
```@repl QI
proj((1/2) * ket(0,2) + (sqrt(3)/2) * ket(1,2))
```
One of the useful operation is matrix vectorization \$\\mathrm{vec}\$, which is a transformation
\$\\mathrm{vec}:\\mathcal{L}(\\mathcal{X,Y})\\to\\mathcal{Y}\\otimes\\mathcal{X}\$ such that
\$\\mathrm{vec}(\\rho)=[a_{1,1},...,a_{m,1}.a_{1,2},...,a_{m,2},...,a_{1,n},....,a_{m,n}]^T\$,
where \$\\rho=[\\rho_{i,j}]\$.
```@repl QI
vec(ketbra(0,1,2))
```
In the opposite to \$\\mathrm{vec}\$ transformation, we have reshaping map \$\\mathrm{res}:\\mathcal{L}(\\mathcal{X,Y})\\to\\mathcal{Y}\\otimes\\mathcal{X}\$, which transform matrix
\$\\rho\$ into a vector row by row. More precisely, for dyadic operators \$|\\psi\\rangle\\langle\\phi|\$, where \$|\\psi\\rangle \\in \\mathcal{X}\$, \$|\\phi\\rangle
\\in \\mathcal{Y}\$ operation \$\\mathrm{vec}\$ is defined as \$\\mathrm{vec}(|\\psi\\rangle\\langle\\phi|)=|\\psi\\rangle|\\overline{\\phi}\\rangle\$ and can be uniquely extend the definition to the whole space \$\\mathcal{L}(\\mathcal{X,Y})\$ by
linearity. It is easy to check that \$\\mathrm{res}(\\rho)=\\mathrm{vec}(\\rho^T)\$.
```@repl QI
res(ketbra(0,1,2))
```

Inverse operation to \$\\mathrm{res}\$ is a \$\\mathrm{unres}:\\mathcal{Y}\\otimes\\mathcal{X}\\to\\mathcal{L}(\\mathcal{X,Y}) \$ map, which transforms
the vector into a matrix. It means that \$\\rho=\\mathrm{unres}(\\mathrm{res}(\\rho))\$.
```@repl QI
unres(res(ketbra(0,1,2)))
```

## Channels
In the most general case evolution of the quantum system can be described
using the notion of a [*quantum channel*](https://en.wikipedia.org/wiki/Quantum_channel).
First, introduce a *superoperator* as a linear mapping acting on linear operators \$\\mathcal{L}(\\mathcal{X})\$
and transforming them into operators acting on \$\\mathcal{L}(\\mathcal{Y})\$. The set of all such mapping wiil be denoted by \$\\mathcal{T}(\\mathcal{X},\\mathcal{Y})\$ and \$\\mathcal{T}(\\mathcal{X}):=\\mathcal{T}(\\mathcal{X},\\mathcal{X})\$. In mathematical terms,
a quantum channel is a superoperator \$\\Phi:\\mathcal{L}(\\mathcal{X})\\to\\mathcal{L}(\\mathcal{Y})\$
that is *trace-preserving* (\$\\forall \\rho\\in\\mathcal{L}(\\mathcal{X})\\quad \\mathrm{Tr}(\\Phi(\\rho))=\\mathrm{Tr}(\\rho)\$)
and *completely positive* (\$\\forall \\mathcal{Z} \\forall \\rho \\in \\mathcal{L}(\\mathcal{X\\otimes Z}), \\rho\\geq 0, \\quad \\Phi\\otimes\\mathbb{I}_{\\mathcal{L}(\\mathcal{X})}(\\rho)\\geq 0\$).
The the product of the given super-operators \$\\Phi_1\\in\\mathcal{T}(\\mathcal{X_1},\\mathcal{Y_1})\$, \$\\Phi_2\\in\\mathcal{T}(\\mathcal{X_2},\\mathcal{Y_2})\$ is a mapping \$\\Phi_1\\otimes\\Phi_2\\in\\mathcal{T}(\\mathcal{X_1}\\otimes\\mathcal{X_2},\\mathcal{Y_1}\\otimes\\mathcal{Y_2})\$ that satisfies \$(\\Phi_1\\otimes\\Phi_2)(\\rho_1\\otimes\\rho_2)=\\Phi_1(\\rho_1)\\otimes\\Phi_2(\\rho_2)\$.

According to Kraus' theorem, any completely positive trace-preserving (CPTP) map \$\\Phi\$ can always be written as
\$\\Phi(\\rho)=\\sum_{i=0}^rK_i \\rho K_i^\\dagger\$ for some set of operators \$\\{K_i\\}_{i=0}^r\$ satisfying \$\\sum_i K_i^\\dagger K_i = \\mathbb{I}\$, where \$r\$ is the rank of super-operator \$\\Phi\$. Another way to represent the quantum channel is based on Choi-Jamiołkowski isomorphism. Consider mapping \$J:\\mathcal{T}(\\mathcal{X,Y})\\to\\mathcal{L}(\\mathcal{Y}\\otimes\\mathcal{X})\$ such that \$J(\\Phi)=(\\Phi\\otimes\\mathbb{I}_{\\mathcal{L}(\\mathcal{X})})(\\mathrm{res}(\\mathbb{I}_{\\mathcal{X}})\\mathrm{res}(\\mathbb{I}_{\\mathcal{X}})^\\dagger)\$. Equivalently \$J(\\Phi)=\\sum_{i,j=0}^{\\mathrm{dim(\\mathcal{X})-1}}\\Phi(|i\\rangle\\langle j|)\\otimes|i\\rangle\\langle j|\$. The action of a superoperator in the Choi representation,also called dynamical matrix of \$\\Phi\$, is given by \$\\Phi(\\rho)=\\mathrm{Tr}_\\mathcal{X}(J(\\Phi)(\\mathbb{I}_\\mathcal{Y}\\otimes\\rho^T))\$. Last representation of quantum channel implemented in `QI.jl` is  Stinespring representation. Supoose that \$A\\in\\mathcal{L}(\\mathcal{X},\\mathcal{Y}\otimes\\mathcal{Z})\$, then \$\\Phi(\\rho)=\\mathrm{Tr}_\\mathcal{Z}(A\\rho A^\\dagger)\$.

### Constructors
Channel objects can be constructed from matrices that represents them. As shown in the following listing
```@repl QI
γ=0.4
K0 = Matrix([1 0; 0 sqrt(1-γ)])
K1 = Matrix([0 sqrt(γ); 0 0])

Φ = KrausOperators([K0,K1])

iscptp(Φ)
istp(Φ)
istni(Φ)
iscptp(Φ)

```
### Conversion
Conversions between all quantum channel types are implemented. The user is not
limited by any single channel representation and can transform between
representations he finds the most efficient or suitable for his purpose.

```@repl QI
Ψ1 = convert(SuperOperator{Matrix{Float64}}, Φ)

Ψ2 = convert(DynamicalMatrix{Matrix{Float64}}, Φ)

Ψ3 = convert(Stinespring{Matrix{Float64}}, Φ)
```

### Application
Channels can act on pure and mixed states as represented by vectors and matrices. Channels
are callable and therefore mimic application of a~function on a~quantum state.

```@repl QI
ψ=(1/sqrt(2)) * (ket(0,2) - ket(1,2))
ρ1=ψ * ψ'

Φ(ρ1)
Ψ1(ρ1)
```

### Composition
Channels can be composed in parallel or in sequence. Composition in parallel is done using
`kron()` function or overloaded \$\\otimes\$ operator. Composition in sequence can
be done in two ways either by using `Julia` build in function composition operator \$(f\\circ g)(\\cdot)=f(g)(\\cdot)\$ or by using multiplication of objects inheriting from `AbstractQuantumOperation{T}` abstract type.

```@repl QI
ϕ=(1/sqrt(2)) * (ket(0,2) - ket(1,2))
ρ2=ϕ * ϕ'

(Φ⊗Φ)(ρ1⊗ρ2)
(Ψ1∘Ψ2)(ρ1)
```

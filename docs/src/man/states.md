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
ket(1,8)
ket(1,8,sparse=true)
```
By \$\\langle\\psi|\$ the
row vector dual to \$|\\psi\\rangle\$ is denoted. It means that \$|\\psi\\rangle=\\langle\\psi|^\\dagger\$, where the symbol \${}^\\dagger\$ denotes the
Hermitian conjugation.
```@repl QI
bra(0,2)  
bra(1,4)
bra(1,4,sparse=true)
```
The form \$|{\\psi}\\rangle\\langle{\\phi}|\$
denotes outer product of \$|{\\psi}\\rangle\$ and \$\\langle{\\phi}|\$ from \$\\mathcal{L}(\\mathcal{X})\$, where  \$\\mathcal{L}(\\mathcal{X})\$ is the set of linear operators
from \$\\mathcal{X}\$ to \$\\mathcal{X}\$. In the case of \$\\mathcal{L}(\\mathcal{X,Y})\$ we have as set of linear operators
from \$\\mathcal{X}\$ to \$\\mathcal{Y}\$ and we assume that \$\\mathcal{L}(\\mathcal{X}) :=\\mathcal{L}(\\mathcal{X,X})\$.
```@repl QI
ketbra(0,1,2)  
ketbra(1,3,4)
ketbra(1,3,4,sparse=true)
```
Especially, \$|{\\psi}\\rangle\\langle{\\psi}|\$ is a rank-one projection operator called as *pure states*. Generally, any [*quantum state*](https://en.wikipedia.org/wiki/Qubit) \$\\rho\$ can be expressed as \$\\rho=\\sum_{i=0}^n q_i |i\\rangle\\langle i|\$, where \$\\sum_{i=0}^n q_i=1\$. Notice that \$\rho\$ is a trace-one positive semi-definite
linear operator *i.e.*: \$\\rho=\\rho^\\dagger\$, \$\\rho\\geq 0\$
and \$\\mathrm{tr}{\\rho}=1\$.
```@repl QI
proj(ket(1,2))
proj(ket(3,4,sparse=true))
```
One of the useful operation is matrix vectorization \$\\mathrm{vec}\$, which is a transformation
\$\\mathrm{vec}:\\mathcal{L}(\\mathcal{X,Y})\\to\\mathcal{Y}\\otimes\\mathcal{X}\$ such that
\$\\mathrm{vec}(\\rho)=[a_{1,1},...,a_{m,1}.a_{1,2},...,a_{m,2},...,a_{1,n},....,a_{m,n}]^T\$,
where \$\\rho=[\\rho_{i,j}]\$.
```@repl QI
vec(ketbra(0,1,2))
```
In the opposite to \$\\mathrm{vec}\$ transformation, we have reshaping map \$\\mathrm{res}:\\mathcal{L}(\\mathcal{X,Y})\\to\\mathcal{X}\\otimes\\mathcal{Y}\$, which transform matrix
\$\\rho\$ into a vector row by row. It is easy to check that \$\\mathrm{res}(\\rho)=\\mathrm{vec}(\\rho^T)\$.
```@repl QI
res(ketbra(0,1,2))
```

Inverse operation to \$\\mathrm{res}\$ is a \$\\mathrm{unres}\$ map, which transforms
the vector into a matrix. It means that \$\\rho=\\mathrm{unres}(\\mathrm{res}(\\rho))\$.
```@repl QI
unres(res(ketbra(0,1,2)))
```

## Channels
In the most general case evolution of the quantum system can be described
using the notion of a [*quantum channel*](https://en.wikipedia.org/wiki/Quantum_channel).
First, introduce a *superoperator* as a linear mapping acting on linear operators \$\\mathcal{L}(\\mathcal{X})\$
and transforming them into operators acting on \$\\mathcal{L}(\\mathcal{Y})\$. In mathematical terms,
a quantum channel is a superoperator \$\\Phi:\\mathcal{L}(\\mathcal{X})\\to\\mathcal{L}(\\mathcal{Y})\$
that is *trace-preserving* (\$\\forall \\rho\\in\\mathcal{L}(\\mathcal{X})\\quad \\mathrm{Tr}(\\Phi(\\rho))=\\mathrm{Tr}(\\rho)\$)
and *completely positive* (\$\\forall \\mathcal{Z} \\forall \\rho \\in \\mathcal{L}(\\mathcal{X\\otimes Z}), \\rho\\geq 0, \\quad \\Phi\\otimes\\mathbb{I}_{\\mathcal{L}(\\mathcal{X})}(\\rho)\\geq 0\$).

```@repl QI
γ=0.4
K0 = Matrix([1 0; 0 sqrt(1-γ)])
K1 = Matrix([0 sqrt(γ); 0 0])
ψ=(1/sqrt(2)) * (ket(0,2) - ket(1,2))
ρ=ψ * ψ'

K = [K0, K1]
Φ = SuperOperator([K0,K1])
println(Φ)
Φ(ρ)

```

According to Kraus' theorem, any completely positive trace-preserving (CPTP) map \$\\Phi\$ can always be written as
\$\\Phi(\\rho)=\\sum_{i}K_i \\rho K_i^\\dagger\$ for some set of operators \$\\{K_i\\}_i\$ satisfying \$\\sum_i K_i^\\dagger K_i = \\mathbb{I}\$.

```@repl QI
Φ = KrausOperators([K0,K1])
println(Φ)

iscptp(Φ)

Φ(ρ)
```

Anotother way to representation of quantum channel is

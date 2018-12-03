```@setup QuantumInformation
using QuantumInformation
```

# Measurement
Measurement is modeled in two ways:
* as Positive Operator Valued Measures (POVMs),
* measurements with post-selection.
In both cases a measurement is treated as a special case of a quantum channel
(operation).

## Positive Operator Valued Measure measurement
A POVM measurement is
defined as follows. Let $\mu:\Gamma\to\mathrm{P}(\mathcal{X})$ be a mapping from
a finite alphabet of measurement outcomes to the set of linear positive
operators. If $\sum_{\xi\in\Gamma} {\mu(\xi)=\mathbb{I}_{\mathcal{X}}}$ then
$\mu$ is a POVM measurement. The set of positive semi-definite linear operators
is defined as $\mathrm{P}(\mathcal{X})=\{X\in \mathrm{L}(\mathcal{X}):
\langle\psi|X|\psi\rangle\geq 0 \text{ for all } |\psi\rangle\in\mathcal{X}\}$. POVM
measurement models the situation where a quantum object is destroyed during the
measurement process and quantum state after the measurement does not exists.

We model POVM measurement as a channel
$\theta:\mathrm{L}(\mathcal{X})\to \mathrm{L}(\mathcal{Y})$, where
$\mathcal{Y}=\mathrm{span}\{|\xi\rangle\}_{\xi\in\Gamma}$ such that $\theta(\rho)
= \sum_{\xi\in\Gamma} \mathrm{Tr}(\rho\, \mu(\xi))|\xi\rangle\langle\xi|$. This channel transforms
the measured quantum state into a classical state (diagonal matrix) containing
probabilities of measuring given outcomes. Note that in **QuantumInformation.jl**
$\Gamma=\{1,2,\ldots,|\Gamma|\}$ and POVM measurements are represented by the
type `POVMMeasurement{T} <: AbstractQuantumOperation{T} where
T<:AbstractMatrix{<:Number}`.
Predicate function `ispovm()` verifies whether a list of matrices is a
proper POVM.

```@repl QuantumInformation
ρ=proj(1.0/sqrt(2)*(ket(1,3)+ket(3,3)))
E0 = proj(ket(1,3))
E1 = proj(ket(2,3))+proj(ket(3,3))

M = POVMMeasurement([E0,E1])

ispovm(M)
M(ρ)
```
## Measurement with post-selection

When a quantum system after being
measured is not destroyed one can be interested in its state after the
measurement. This state depends on the measurement outcome. In this case the
measurement process is defined in the following way.

Let $\mu:\Gamma\to \mathrm{L}(\mathcal{X}, \mathcal{Y})$ be a mapping from a finite
set of measurement outcomes to set of linear operators called effects. If
$\sum_{\xi\in\Gamma} {\mu(\xi)^\dagger \mu(\xi)=\mathbb{I}_{\mathcal{X}}}$ then
$\mu$ is a quantum measurement. Given outcome $\xi$ was obtained, the state
before the measurement, $\rho$, is transformed into sub-normalized quantum state
$\rho_\xi=\mu(\xi)\rho\mu(\xi)^\dagger$. The outcome $\xi$ will be obtained with
probability $\mathrm{Tr}(\rho_\xi)$.
```@repl QuantumInformation
PM = PostSelectionMeasurement(E1)
iseffect(PM)
PM(ρ)
```
In **QuantumInformation** this kind of measurement is modeled as CP-TNI map with single Kraus operator $\mu(\xi)$ and represented as
`PostSelectionMeasurement{T} <: AbstractQuantumOperation{T} where T<:AbstractMatrix{<:Number}`. Measurement types can be composed and converted to Kraus operators,
superoperators, Stinespring representation operators, and dynamical matrices.
```@repl QuantumInformation
α = 0.3
K0 = ComplexF64[0 0 sqrt(α); 0 1 0; 0 0 0]
K1 = ComplexF64[1 0 0; 0 0 0; 0 0 sqrt(1 - α)]
Φ = KrausOperators([K0,K1])

ρ=proj(1.0/sqrt(2)*(ket(1,3)+ket(3,3)))

(PM∘Φ)(ρ)
```

```@setup QuantumInformation
using QuantumInformation
using LinearAlgebra
```

# Channels

This section describes functions for working with Quantum Channels, including application, representation, and property checking.

## Applying Channels

Quantum channels can be applied to quantum states (density matrices) and state vectors (kets). The `applychannel` function handles various channel representations such as `KrausOperators`, `SuperOperator`, `CholJamiolkowskiMatrices` (as `DynamicalMatrix`), `Stinespring`, `IdentityChannel`, and `UnitaryChannel`.

```@repl QuantumInformation
ρ = [0.5 0; 0 0.5]
Φ = IdentityChannel(2)
applychannel(Φ, ρ)

U = UnitaryChannel([0 1; 1 0])
applychannel(U, ρ)
```

Also, channels object are callable
```@repl QuantumInformation
Φ(ρ)
```

## Channel Representations

Channels can be represented in various bases. The `channelbasis` function generates a basis for quantum channels, which can then be used with `represent` to find the coefficients of a channel in that basis, or `combine` to reconstruct the channel from coefficients.

```@repl QuantumInformation
# Create a channel basis for qubits (idim=2, odim=2)
basis = channelbasis(2)

# Represent a channel (e.g., Identity) in this basis
id_op = vec(Matrix{ComplexF64}(I, 4, 4))
coeffs = represent(basis, reshape(id_op, 4, 4)) # Example for illustration

# Reconstruct
# op = combine(basis, coeffs)
```

## Predicates

We provide a set of predicates to check properties of quantum operations.

### CP (Completely Positive)
- `iscp(Φ)`: Checks if `Φ` is Completely Positive.

### TP (Trace Preserving) and TNI (Trace Non-Increasing)
- `istp(Φ)`: Checks if `Φ` is Trace Preserving.
- `istni(Φ)`: Checks if `Φ` is Trace Non-Increasing.

### CPTP (Quantum Channel)
- `iscptp(Φ)`: Checks if `Φ` is both CP and TP.
- `iscptni(Φ)`: Checks if `Φ` is both CP and TNI.

### Measurements
- `ispovm(Φ)`: Checks if a set of operators forms a POVM.
- `iseffect(Φ)`: Checks if an operator is a valid quantum effect.

```@repl QuantumInformation
K = KrausOperators([[1 0; 0 1]])
iscptp(K)
```

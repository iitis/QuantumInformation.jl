
```@meta
Author = "Piotr Gawron, Dariusz Kurzyk, Łukasz Pawela"
```

# Home

A Julia package for numerical computation in quantum information theory.

Numerical investigations are prevalent in quantum information theory. Numerical experiments can be used to find counter examples for theorems, to test hypotheses or to gain insight about quantum objects and operations.

Our goal while designing **QI.jl** library was to follow principles presented in book "Geometry of Quantum States'' [[1]](@ref refs). We work with column vectors reprinting kets and row vectors representing bras. We fix our basis to the computational one. Density matrices and quantum channels are represented as two dimensional arrays in the same fixed basis. This approach allows us to obtain low level complexity of our code, high flexibility and good computational efficiency. The design choices where highly motivated by the properties of the language in which the our library was implemented, namely
[Julia](https://julialang.org/) [[2]](@ref refs).

## Package features


## Credits


## [References](@id refs)

[1] I. Bengtsson, K. Życzkowski, *Geometry of Quantum States: An Introduction to Quantum Entanglement*, Cambridge University Press (2008).

[2] J. Bezanson, S. Karpinski, V. B. Shah, A. Edelman, *Julia: A fast dynamic language for technical computing*,
[preprint](https://arxiv.org/pdf/1209.5145.pdf).

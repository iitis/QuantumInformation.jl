var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "Author = \"Piotr Gawron, Dariusz Kurzyk, Łukasz Pawela\""
},

{
    "location": "#Home-1",
    "page": "Home",
    "title": "Home",
    "category": "section",
    "text": "A Julia package for numerical computation in quantum information theory.Numerical investigations are prevalent in quantum information theory. Numerical experiments can be used to find counter examples for theorems, to test hypotheses or to gain insight about quantum objects and operations.Our goal while designing QuantumInformation.jl library was to follow principles presented in book \"Geometry of Quantum States\'\' [1]. We work with column vectors reprinting kets and row vectors representing bras. We fix our basis to the computational one. Density matrices and quantum channels are represented as two dimensional arrays in the same fixed basis. This approach allows us to obtain low level complexity of our code, high flexibility and good computational efficiency. The design choices where highly motivated by the properties of the language in which the our library was implemented, namely Julia [2]."
},

{
    "location": "#Package-features-1",
    "page": "Home",
    "title": "Package features",
    "category": "section",
    "text": "The purpose of QuantumInformation.jl library is to provide functions to:creating and analyzing quantumstates,manipulating them with quantum channels\ncalculating functionals on these objects, i.e. trace norm, diamond norm, entropy, fidelity,\napplication of random matrix theory in quantuminformation processing."
},

{
    "location": "#refs-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "[1] I. Bengtsson, K. Życzkowski, Geometry of Quantum States: An Introduction to Quantum Entanglement, Cambridge University Press (2008).[2] J. Bezanson, S. Karpinski, V. B. Shah, A. Edelman, Julia: A fast dynamic language for technical computing, preprint."
},

{
    "location": "man/quickstart/#",
    "page": "Quickstart",
    "title": "Quickstart",
    "category": "page",
    "text": ""
},

{
    "location": "man/quickstart/#manual-1",
    "page": "Quickstart",
    "title": "Quickstart",
    "category": "section",
    "text": "If this is your first use of the Julia Language, please read the Documentation. The latest release of QuantumInformation.jl can be installed from the Julia REPL prompt withjulia> Pkg.add(\"QuantumInformation.jl\")Package QuantumInformation can be loaded viajulia> using QuantumInformationThis manual provides a introduction explaining how use this package for numerical computations in quantum information theory."
},

{
    "location": "man/vectors/#",
    "page": "Vectors and matrices in Julia",
    "title": "Vectors and matrices in Julia",
    "category": "page",
    "text": ""
},

{
    "location": "man/vectors/#Vectors-and-matrices-in-Julia-1",
    "page": "Vectors and matrices in Julia",
    "title": "Vectors and matrices in Julia",
    "category": "section",
    "text": "<!– # Vectors and matrices in JuliaA basic construction of vector in Julia creates a full one-index array containing elements of a number type as presented below.julia> x = [0.0, 1.0im]\n2-element Array{Complex{Float64},1}:\n 0.0+0.0im\n 0.0+1.0imA transposition of a column vector return an object of type RowVector as shown belowjulia> xt = transpose(x)\n1×2 RowVector{Complex{Float64},Array{Complex{Float64},1}}:\n 0.0+0.0im  0.0+1.0imWhile a hermitian conjugate of the same vector returns a RowVector parametrized by the type ConjArrayjulia> xc = [0.0, 1.0im]\'\n1×2 RowVector{Complex{Float64},ConjArray{Complex{Float64},1,Array{Complex{Float64},1}}}:\n 0.0-0.0im  0.0-1.0imValues of variables xt and xc are views of the value of variable x. When used the column and row vectors behave like bra and kets, for example xc*x denotes the inner product of bra xc and ket x, while x*xc denotes its outer product resulting in a two-index array.The linear algebra library in Julia provides standard operations on matrices and vectors that are designed to take in to the account the types of objects. –>"
},

{
    "location": "man/states/#",
    "page": "States",
    "title": "States",
    "category": "page",
    "text": ""
},

{
    "location": "man/states/#States-1",
    "page": "States",
    "title": "States",
    "category": "section",
    "text": "<!– ```@setup QuantumInformation Pkg.add(\"QuantumInformation\") importall QuantumInformation\n# States and channels\nIn this and the following sections we will denote complex Euclidean spaces\n$\\mathbb{C}^d$ with $\\mathcal{X}$, $\\mathcal{Y}$, $\\mathcal{Z}$ etc. When needed the dimension of a space $\\mathcal{X}$\nwill be denoted $\\mathrm{dim}(\\mathcal{X})$. The set of matrices transforming vectors\nfrom $\\mathcal{X}$ to $\\mathcal{Y}$ will be denoted $\\mathrm{L}(\\mathcal{X}, \\mathcal{Y})$. For simplicity we will\nwrite $\\mathrm{L}(\\mathcal{X}) \\equiv \\mathrm{L}(\\mathcal{X}, \\mathcal{X})$.\n\n\n## States\nBy\n\\$|\\\\psi\\\\rangle\\in\\\\mathcal{X}\\$ we denote a normed column vector. Notice that\nany \\$|\\\\psi\\\\rangle\\$ can be expressed as \\$|\\\\psi\\\\rangle=\\\\sum_{i=0}^n \\\\alpha_i |i\\\\rangle\\$, where \\$\\\\sum_{i=0}^n |\\\\alpha_i|^2=1\\$ and set \\$\\\\{|i\\\\rangle\\\\}_{i=0}^n\\$ is a computational basis.@repl QuantumInformation ket(0,2)   (1/sqrt(2)) * (ket(0,2) + ket(1,2))By \\$\\\\langle\\\\psi|\\$ the\nrow vector dual to \\$|\\\\psi\\\\rangle\\$ is denoted. It means that \\$|\\\\psi\\\\rangle=\\\\langle\\\\psi|^\\\\dagger\\$, where the symbol \\${}^\\\\dagger\\$ denotes the\nHermitian conjugation.@repl QuantumInformation bra(1,3)  \nThrough the \\$\\\\langle\\\\psi|\\\\phi\\\\rangle\\$ is denoted a inner product and the\nnorm defined as \\$\\\\|\\\\phi\\\\rangle\\\\|=\\\\sqrt{\\\\langle\\\\psi|\\\\psi\\\\rangle}\\$.@repl QuantumInformation ψ=(1/sqrt(2)) * (ket(0,2) + ket(1,2)) ϕ=(1/2) * ket(0,2) + (sqrt(3)/2) * ket(1,2)ϕ\'*ψsqrt(ϕ\'*ϕ)\nThe form \\$|{\\\\psi}\\\\rangle\\\\langle{\\\\phi}|\\$\ndenotes outer product of \\$|{\\\\psi}\\\\rangle\\$ and \\$\\\\langle{\\\\phi}|\\$ from \\$\\\\mathrm{L}(\\\\mathcal{X})\\$.@repl QuantumInformation ketbra(1,3,4)Specifically, \\$|{\\\\psi}\\\\rangle\\\\langle{\\\\psi}|\\$ is a rank-one projection operator called as *pure states*. Generally, any [*quantum state*](https://en.wikipedia.org/wiki/Qubit) \\$\\\\rho\\$ can be expressed as \\$\\\\rho=\\\\sum_{i=0}^n q_i |i\\\\rangle\\\\langle i|\\$, where \\$\\\\sum_{i=0}^n q_i=1\\$. Notice that \\$\\rho\\$ is a trace-one positive semi-definite\nlinear operator *i.e.*: \\$\\\\rho=\\\\rho^\\\\dagger\\$, \\$\\\\rho\\\\geq 0\\$\nand \\$\\\\mathrm{tr}{\\\\rho}=1\\$.@repl QuantumInformation proj(ψ)In the opposite to \\$\\\\mathrm{vec}\\$ transformation, we have reshaping map \\$\\\\mathrm{res}:\\\\mathrm{L}(\\\\mathcal{X,Y})\\\\to\\\\mathcal{Y}\\\\otimes\\\\mathcal{X}\\$, which transform matrix\n\\$\\\\rho\\$ into a vector row by row. More precisely, for dyadic operators \\$|\\\\psi\\\\rangle\\\\langle\\\\phi|\\$, where \\$|\\\\psi\\\\rangle \\\\in \\\\mathcal{X}\\$, \\$|\\\\phi\\\\rangle\n\\\\in \\\\mathcal{Y}\\$ operation \\$\\\\mathrm{vec}\\$ is defined as \\$\\\\mathrm{vec}(|\\\\psi\\\\rangle\\\\langle\\\\phi|)=|\\\\psi\\\\rangle|\\\\overline{\\\\phi}\\\\rangle\\$ and can be uniquely extend the definition to the whole space \\$\\\\mathrm{L}(\\\\mathcal{X,Y})\\$ by\nlinearity. It is easy to check that \\$\\\\mathrm{res}(\\\\rho)=\\\\mathrm{vec}(\\\\rho^T)\\$.@repl QuantumInformation res(ketbra(0,1,2))\nInverse operation to \\$\\\\mathrm{res}\\$ is a \\$\\\\mathrm{unres}:\\\\mathcal{Y}\\\\otimes\\\\mathcal{X}\\\\to \\\\mathrm{L}(\\\\mathcal{X,Y}) \\$ map, which transforms\nthe vector into a matrix. It means that \\$\\\\rho=\\\\mathrm{unres}(\\\\mathrm{res}(\\\\rho))\\$.@repl QuantumInformation unres(res(ketbra(0,1,2)))Let us recall that trace is a mapping \\$\\\\mathrm{Tr}:\\\\mathrm{L}(\\\\mathcal{X})\\\\ \\\\to \\\\mathbb{C},\\$ given by \\$\\\\mathrm{Tr}:\\\\rho\\\\mapsto\\\\sum_{i=1}^{\\\\mathrm{dim}(\\\\mathcal{X})}\\\\langle e_i|\\\\rho|e_i\\\\rangle\\$, where \\$\\\\{|e_i\\\\rangle \\\\}\\$ is an orthonormal basis of \\$\\\\mathcal{X}\\$. According this, [*partial trace*](https://en.wikipedia.org/wiki/Partial_trace) is a mapping \\$\\\\mathrm{Tr}_{\\\\mathcal{X}}: \\\\mathrm{L}(\\\\mathcal{X}\\\\otimes\\\\mathcal{Y}) \\\\to \\\\mathrm{L}(\\\\mathcal{Y})\\$ such that \\$\\\\mathrm{Tr}_{\\\\mathcal{X}}: \\\\rho_A\\\\otimes \\\\rho_B \\\\mapsto \\\\rho_B \\\\mathrm{Tr}(\\\\rho_A)\\$, where \\$\\\\rho_A\\\\in \\\\mathrm{L}(\\\\mathcal{X})\\$, \\$\\\\rho_B\\\\in \\\\mathrm{L}(\\\\mathcal{Y})\\$.\nAs this is a linear map, it may be uniquely extended to the case of operators\nwhich are not in a tensor product form.@repl QuantumInformation ρ = [0.25 0.25im; -0.25im 0.75] σ = [0.4 0.1im; -0.1im 0.6]ptrace(ρ ⊗ σ, [2, 2], [2,])\nMatrix transposition is a mapping \\${}^T:\\\\mathrm{L}(\\\\mathcal{X,Y}) \\\\to \\\\mathrm{L}(\\\\mathcal{Y,X})\\$ such that \\$\\\\rho_{ij}^T = \\\\rho_{ji}\\$, where \\$\\\\rho_{ij}\\$ is a \\$i\\$-th row, \\$j\\$-th column element of matrix \\$\\\\rho\\$. Following this, we may introduce *partial transposition* \\${}^{T_B}:\\\\mathrm{L}(\\\\mathcal{X_A,Y_A}\\\\otimes\\\\mathcal{X_B,Y_B}) \\\\to \\\\mathrm{L}(\\\\mathcal{Y_A,X_A}\\\\otimes\\\\mathcal{X_B,Y_B})\\$,\nwhich for product state \\$\\\\rho_A\\\\otimes\\\\rho_B\\$ is given by \\${}^{T_A}: \\\\rho_A\\\\otimes\\\\rho_B\\\\mapsto\\\\rho_A^T\\\\otimes\\\\rho_B\\$.\nThe definition of partial trace can be extended for all operators by trace linearity.@repl QuantumInformation ptranspose(ρ ⊗ σ, [2, 2], [1])\nFor given multiindexed matrix \\$\\\\rho_{(m,\\\\mu),(n,\\\\nu)}=\\\\langle m|\\\\langle\\\\mu|\\\\rho|n\\\\rangle\\\\nu\\\\rangle\\$,\nreshuffle operation is defined as \\$\\\\rho^R_{(m,\\\\mu),(n,\\\\nu)}=\\\\rho_{(m,n),(\\\\mu,\\\\nu)}\\$.@repl QuantumInformation reshuffle(ρ ⊗ σ)\n## Channels\nIn the most general case evolution of the quantum system can be described\nusing the notion of a [*quantum channel*](https://en.wikipedia.org/wiki/Quantum_channel).\nFirst, introduce a *superoperator* as a linear mapping acting on linear operators \\$\\\\mathrm{L}(\\\\mathcal{X})\\$\nand transforming them into operators acting on \\$\\\\mathrm{L}(\\\\mathcal{Y})\\$. The set of all such mapping will be denoted by \\$\\\\mathrm{T}(\\\\mathcal{X},\\\\mathcal{Y})\\$ and \\$\\\\mathrm{T}(\\\\mathcal{X}):=\\\\mathrm{T}(\\\\mathcal{X},\\\\mathcal{X})\\$. In mathematical terms,\na quantum channel is a superoperator \\$\\\\Phi:\\\\mathrm{L}(\\\\mathcal{X})\\\\to \\\\mathrm{L}(\\\\mathcal{Y})\\$\nthat is *trace-preserving* (\\$\\\\forall \\\\rho\\\\in \\\\mathrm{L}(\\\\mathcal{X})\\\\quad \\\\mathrm{Tr}(\\\\Phi(\\\\rho))=\\\\mathrm{Tr}(\\\\rho)\\$)\nand *completely positive* (\\$\\\\forall \\\\mathcal{Z} \\\\forall \\\\rho \\\\in \\\\mathrm{L}(\\\\mathcal{X\\\\otimes Z}), \\\\rho\\\\geq 0, \\\\quad \\\\Phi\\\\otimes\\\\mathbb{I}_{\\\\mathrm{L}(\\\\mathcal{X})}(\\\\rho)\\\\geq 0\\$).\nThe the product of the given super-operators \\$\\\\Phi_1\\\\in \\\\mathrm{T}(\\\\mathcal{X_1},\\\\mathcal{Y_1})\\$, \\$\\\\Phi_2\\\\in \\\\mathrm{T}(\\\\mathcal{X_2},\\\\mathcal{Y_2})\\$ is a mapping \\$\\\\Phi_1\\\\otimes\\\\Phi_2\\\\in T(\\\\mathcal{X}_1\\\\otimes\\\\mathcal{X}_2,\\\\mathcal{Y}_1\\\\otimes\\\\mathcal{Y}_2)\\$ that satisfies \\$(\\\\Phi_1\\\\otimes\\\\Phi_2)(\\\\rho_1\\\\otimes\\\\rho_2)=\\\\Phi_1(\\\\rho_1)\\\\otimes\\\\Phi_2(\\\\rho_2)\\$.\n\nAccording to Kraus\' theorem, any completely positive trace-preserving (CPTP) map \\$\\\\Phi\\$ can always be written as\n\\$\\\\Phi(\\\\rho)=\\\\sum_{i=0}^rK_i \\\\rho K_i^\\\\dagger\\$ for some set of operators \\$\\\\{K_i\\\\}_{i=0}^r\\$ satisfying \\$\\\\sum_i K_i^\\\\dagger K_i = \\\\mathbb{I}\\$, where \\$r\\$ is the rank of super-operator \\$\\\\Phi\\$. Another way to represent the quantum channel is based on Choi-Jamiołkowski isomorphism. Consider mapping \\$J:\\\\mathrm{T}(\\\\mathcal{X,Y})\\\\to \\\\mathrm{L}(\\\\mathcal{Y}\\\\otimes\\\\mathcal{X})\\$ such that \\$J(\\\\Phi)=(\\\\Phi\\\\otimes\\\\mathbb{I}_{\\\\mathrm{L}(\\\\mathcal{X})})(\\\\mathrm{res}(\\\\mathbb{I}_{\\\\mathcal{X}})\\\\mathrm{res}(\\\\mathbb{I}_{\\\\mathcal{X}})^\\\\dagger)\\$. Equivalently \\$J(\\\\Phi)=\\\\sum_{i,j=0}^{\\\\mathrm{dim(\\\\mathcal{X})-1}}\\\\Phi(|i\\\\rangle\\\\langle j|)\\\\otimes|i\\\\rangle\\\\langle j|\\$. The action of a superoperator in the Choi representation,also called dynamical matrix of \\$\\\\Phi\\$, is given by \\$\\\\Phi(\\\\rho)=\\\\mathrm{Tr}_\\\\mathcal{X}(J(\\\\Phi)(\\\\mathbb{I}_\\\\mathcal{Y}\\\\otimes\\\\rho^T))\\$. Last representation of quantum channel implemented in `QuantumInformation.jl` is  Stinespring representation. Supoose that \\$A\\\\in \\\\mathrm{L}(\\\\mathcal{X},\\\\mathcal{Y}\\otimes\\\\mathcal{Z})\\$, then \\$\\\\Phi(\\\\rho)=\\\\mathrm{Tr}_\\\\mathcal{Z}(A\\\\rho A^\\\\dagger)\\$.\n\n### Constructors\nChannel objects can be constructed from matrices that represents them. As shown in the following listing@repl QuantumInformation γ=0.4 K0 = Matrix([1 0; 0 sqrt(1-γ)]) K1 = Matrix([0 sqrt(γ); 0 0])Φ = KrausOperators([K0,K1])iscptp(Φ)### Conversion\nConversions between all quantum channel types are implemented. The user is not\nlimited by any single channel representation and can transform between\nrepresentations he finds the most efficient or suitable for his purpose.\n@repl QuantumInformation Ψ1 = convert(SuperOperator{Matrix{ComplexF64}}, Φ)"
},

{
    "location": "man/states/#or-1",
    "page": "States",
    "title": "or",
    "category": "section",
    "text": "SuperOperator{Matrix{ComplexF64}}(Φ)Ψ2 = convert(DynamicalMatrix{Matrix{Float64}}, Φ)#orDynamicalMatrix{Matrix{Float64}}(Φ)Ψ3 = convert(Stinespring{Matrix{Float64}}, Φ)#orStinespring{Matrix{Float64}}(Φ)\n### Application\nChannels can act on pure and mixed states as represented by vectors and matrices. Channels\nare callable and therefore mimic application of a~function on a~quantum state.\n@repl QuantumInformation ρ1=ψ * ψ\'Φ(ρ1) Ψ1(ρ1) Φ(ψ)\n### Composition\nChannels can be composed in parallel or in sequence. Composition in parallel is done using\n`kron()` function or overloaded \\$\\\\otimes\\$ operator. Composition in sequence can\nbe done in two ways either by using `Julia` build in function composition operator \\$(f\\\\circ g)(\\\\cdot)=f(g)(\\\\cdot)\\$ or by using multiplication of objects inheriting from `AbstractQuantumOperation{T}` abstract type.\n@repl QuantumInformation ρ2=ϕ * ϕ\'(Φ⊗Φ)(ρ1⊗ρ2) (Ψ1∘Ψ2)(ρ1)\n## Measurement\nMeasurement is modeled in two ways:\n* as Positive Operator Valued Measure (POVM) based,\n* measurement with post-selection.\nIn both cases a~measurement is treated as a special case of quantum channel\n(operation) respectively as defined below.\n\n### Positive Operator Valued Measure measurement\nA POVM measurement is defined as follows, let a\nmapping from a finite alphabet of measurement outcomes to the set of linear\npositive operators \\$\\\\mu:\\\\Gamma\\\\rightarrow\\\\mathcal{P}(\\\\mathcal{X})\\$\n be given, if \\$\\\\sum_{\\\\xi\\\\in\\\\Gamma} {\\\\mu(\\\\xi)=\\\\mathbb{I}_{\\\\mathcal{X}}}\\$ then \\$\\\\mu\\$ is a POVM measurement. POVM\nmeasurement models the situation where a quantum object is destroyed during the\nmeasurement process and quantum state after the measurement does not exists.\n\nWe model POVM measurement as a channel\n\\$\\\\theta:\\\\mathcal{T}(\\\\mathcal{X})\\\\rightarrow \\\\mathcal{T}(\\\\mathcal{Y})\\$, where\n\\$\\\\mathcal{Y}=\\\\mathrm{span}\\\\{|\\\\xi\\\\rangle\\\\}_{\\\\xi\\\\in\\\\Gamma}\\$ such that \\$\\\\theta(\\\\rho)\n= \\\\sum_{\\\\xi\\\\in\\\\Gamma} \\\\mathrm{tr}(\\\\rho\\\\, \\\\mu(\\\\xi))\\\\|\\\\xi\\\\rangle\\\\langle \\\\xi|\\$. This channel transforms\nthe measured quantum state into a~classical state (diagonal matrix) containing\nprobabilities of reassuring given outcomes. Note that in `QuantumInformation`\n\\$\\\\Gamma=\\\\{1,2,\\\\ldots,|\\\\Gamma|\\\\}\\$ and POVM measurements are represented by the\ntype `POVMMeasurement{T} <: AbstractQuantumOperation{T}` where `T<:AbstractMatrix{<:Number}`.\nPredicate function `ispovm()` verifies whether a~list of matrices is a proper POVM.@repl QuantumInformation ρ=proj(1./sqrt(2)*(ket(0,3)+ket(2,3)))E0 = proj(ket(0,3)) E1 = proj(ket(1,3))+proj(ket(2,3))M = POVMMeasurement([E0,E1])ispovm(M) M(ρ)When a quantum system after being measured is not destroyed one can be\ninterested by its state after the measurement. This state depends on the\nmeasurement outcome. In this case the measurement process is defined in the following way.\n\nLet $\\mu:\\Gamma\\rightarrow\\mathcal{L}(\\mathcal{X}, \\mathcal{Y})$\nbe a mapping from a finite set of measurement outcomes to set of linear operators called effects, then\nif $\\sum_{\\xi\\in\\Gamma} {\\mu(\\xi)^\\dagger \\mu(\\xi)=\\mathbb{I}_{\\mathcal{X}}}$\nthen $\\mu$ is a quantum measurement. Given outcome $\\xi$ was obtained, the state before the measurement $\\rho$\nis transformed into sub-normalized quantum state $\\rho_\\xi=\\mu(\\xi)\\rho\\mu(\\xi)^\\dagger$. The outcome $\\xi$ will be obtained\nwith probability $\\mathrm{Tr}(\\rho_\\xi)$.@repl QuantumInformation PM = PostSelectionMeasurement(E1) iseffect(PM) PM(ρ)In `QuantumInformation` this kind of measurement is modeled as CP-TNI map with single Kraus operator $\\mu(\\xi)$ and represented as\n`PostSelectionMeasurement{T} <: AbstractQuantumOperation{T}` where `T<:AbstractMatrix{<:Number}`. Measurement types can be composed and converted to Kraus operators,\nsuperoperators, Stinespring representation operators, and dynamical matrices.@repl QuantumInformation α = 0.3 K0 = ComplexF64[0 0 sqrt(α); 0 1 0; 0 0 0] K1 = ComplexF64[1 0 0; 0 0 0; 0 0 sqrt(1 - α)] Φ = KrausOperators([K0,K1])ρ=proj(1./sqrt(2)*(ket(0,3)+ket(2,3)))(PM∘Φ)(ρ) ``` –>"
},

{
    "location": "man/functionals/#",
    "page": "Functionals",
    "title": "Functionals",
    "category": "page",
    "text": ""
},

{
    "location": "man/functionals/#Functionals-1",
    "page": "Functionals",
    "title": "Functionals",
    "category": "section",
    "text": "<!– ```@setup QuantumInformation Pkg.add(\"QuantumInformation\") importall QuantumInformation\n# Functionals\n\nLet us \\$\\\\rho, \\\\sigma \\\\in \\\\mathrm{L}(\\\\mathcal{X})\\$. [*Trace norm*](https://www.quantiki.org/wiki/trace-norm) is defined as \\$\\\\|\\\\rho\\\\|_1 = \\\\mathrm{Tr} \\\\sqrt{\\\\rho\\\\rho^\\\\dagger}\\$ and trace distance is defined based on the trace norm as \\$D_1(\\\\rho,\\\\sigma)=\\\\frac{1}{2}\\\\|\\\\rho-\\\\sigma\\\\|_1$.\n@repl QuantumInformation ψ=(1/sqrt(2)) * (ket(0,2) + ket(1,2)) ϕ=(1/2) * ket(0,2) + (sqrt(3)/2) * ket(1,2)ρ=proj(ψ) σ=proj(ϕ)normtrace(ρ) tracedistance(ρ, σ)Now introduce [*Hilbert–Schmidt norm*](https://en.wikipedia.org/wiki/Hilbert%E2%80%93Schmidt_operator) and distance by \\$\\\\|\\\\rho\\\\|_{HS}=\\\\mathrm{Tr}\\\\rho^\\\\dagger \\\\rho\\$ and \\$D_{HS}(\\\\rho,\\\\sigma)=\\\\frac{1}{2}\\\\|\\\\rho-\\\\sigma\\\\|_{HS}$, respectively.@repl QuantumInformation normhs(ρ) hsdistance(ρ, σ)\n[*Fidelity*](https://en.wikipedia.org/wiki/Fidelity_of_quantum_states) is a measure of distance of quantum states. It is an example of a\ndistance measure which is not a metric on the space of quantum states. The\nfidelity of two quantum states \\$\\\\rho, \\\\sigma \\in \\\\mathrm{L}(\\\\mathcal{X})\\$ is given by\n\\$F(\\\\rho,\\\\sigma)=\\\\|\\\\sqrt{\\\\rho}\\\\sqrt{\\\\sigma}\\\\|_1\\$@repl QuantumInformation fidelity_sqrt(ρ, σ) fidelity(ρ, σ) fidelity(ψ, σ) fidelity(ρ, ϕ) fidelity(ψ, ϕ)\n[*Superfidelity*](https://www.quantiki.org/wiki/superfidelity) is an upper bound on the fidelity of two quantum states\nIt is defined by\n\\$G(\\\\rho, \\\\sigma) = \\mathrm{Tr}\\\\rho \\\\sigma + \\\\sqrt{1 - \\mathrm{Tr}\\\\rho^2} \\\\sqrt{1-\\mathrm{Tr} \\\\sigma^2}\\$.\n@repl QuantumInformation superfidelity(ρ, σ)\nIn order to introduce the diamond norm, we first introduce the notion of the\ninduced trace norm. Given \\$\\\\Phi \\\\in \\\\mathrm{T}(\\\\mathcal{X}, \\\\mathcal{Y})\\$ we define its induced trace\nnorm as \\$\\\\| \\\\Phi \\\\|_1 = \\\\mathrm{max} \\\\left\\\\{ \\\\| \\\\Phi(X) \\\\|_1: X \\\\in L(\\\\mathcal{X}), \\\\| X \\\\|_1 \\\\leq 1\n\\\\right\\\\}\\$.\nThe diamond norm of \\$\\\\Phi\\$ is defined as\n\\$\n\\\\| \\\\Phi \\\\|_\\\\diamond = \\\\| \\\\Phi \\\\otimes \\\\mathbb{I} \\\\|_1\n\\$\nOne important property of the diamond norm is that for Hermiticity-preserving\n\\$\\\\Phi \\\\in \\\\mathrm{T}(\\\\mathcal{X}, \\\\mathcal{Y})\\$ we obtain\n\\$\n\\\\| \\\\Phi \\\\|_\\\\diamond = \\\\max \\\\left\\\\{ \\\\left\\\\| (\\\\Phi \\\\otimes \\\\mathbb{I})\n\\\\left(|\\\\psi\\\\rangle\\\\langle\\\\psi| \\\\right )\\\\right\\\\|_1: |\\\\psi\\\\rangle \\\\in \\\\mathcal{X} \\\\otimes \\\\mathcal{Y},\n\\\\langle\\\\psi|\\\\psi\\\\rangle=1 \\\\right\\\\}\\$.\n@repl QuantumInformation K0 = Matrix([1 0; 0 sqrt(1-γ)]) K1 = Matrix([0 sqrt(γ); 0 0])Φ = KrausOperators([K0,K1])L0 = Matrix([1 0; 0 sqrt(1-γ)]) L1 = Matrix([0 0; 0 sqrt(γ)])Ψ = KrausOperators([K0,K1])normdiamond(Φ) diamonddistance(Φ, Ψ)\n[*Shannon entropy*](https://en.wikipedia.org/wiki/Entropy_(information_theory)) is defined as \\$H(\\\\mathrm{p})=-\\\\sum_{i=1}^n p_i\\\\log_2 p_i\\$, where\n\\$\\\\mathrm{p}=[p_1,...,p_n]\\$ is a vector of probabilities.\n@repl QuantumInformation p = vec([0.3 0.2 05])shannonentropy(p) shannonentropy(0.5)\nFor a quantum system described by a state \\$\\\\rho\\$, the [*von Neumann entropy*](https://en.wikipedia.org/wiki/Von_Neumann_entropy) is \\$S(\\\\rho)=-\\\\mathrm{tr} \\\\rho \\\\log \\\\rho\\$.\nLet \\$\\\\lambda_i\\$,  \\$0\\\\geq i < n\\$ be eigenvalues of \\$\\\\rho\\$, then \\$S(\\\\rho)\\$ can be written as \\$S(\\\\rho)=-\\\\sum_{i=0}^{n-1} \\\\lambda_i \\\\log \\\\lambda_i\\$.@repl QuantumInformation ρ = [0.25 0.25im; -0.25im 0.75] σ = [0.4 0.1im; -0.1im 0.6]quantum_entropy(0.4 * ρ + 0.6 * σ)\nOne of the measure of distinguishability between two quantum states is a [*qauntum relative entropy*](https://en.wikipedia.org/wiki/Quantum_relative_entropy), called also Kullback–Leibler divergence, defined as\n\\$S(\\\\rho\\\\|\\\\sigma)=-\\\\mathrm{Tr}\\\\rho\\\\log\\\\sigma + \\\\mathrm{Tr}\\\\rho\\\\log\\\\rho\\$@repl QuantumInformation relativeentropy(ρ, σ) kldivergence(ρ, σ)\nAnother type of measure of distinguishability between two quantum state is [*quantum Jensen–Shannon divergence*](https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence#Quantum_Jensen%E2%80%93Shannon_divergence) given by\n\\$QJS(\\\\rho,\\\\sigma)=S\\\\left(\\\\frac{1}{2}\\\\rho+\\\\frac{1}{2}\\\\sigma\\\\right)-\\\\left(\\\\frac{1}{2}S(\\\\rho)+\\\\frac{1}{2}S(\\\\sigma)\\\\right)\\$.@repl QuantumInformation js_divergence(ρ, σ)\n[*The Bures distance*](https://en.wikipedia.org/wiki/Bures_metric) defines an infinitesimal distance between quantum states, and it is defined as \\$D_B=\\\\sqrt{2(1-\\\\sqrt{F(\\\\rho,\\\\sigma)})}\\$. The value related with Bures distance is Bures angle \\$D_A(\\\\rho,\\\\sigma)=\\\\arccos(\\\\sqrt{F(\\\\rho,\\\\sigma)})\\$@repl QuantumInformation buresdistance(ρ, σ) buresangle(ρ, σ)\nOne of the entanglement measure is [*negativity*](https://en.wikipedia.org/wiki/Negativity_(quantum_mechanics)) defined as \\$\\\\mathcal{N}(\\\\rho)=\\\\frac{\\\\|\\\\rho^{T_A}\\\\|_1-1}{2}\\$.\ndistance is Bures angle \\$D_A(\\\\rho,\\\\sigma)=\\\\arccos(\\\\sqrt{F(\\\\rho,\\\\sigma)})\\$@repl QuantumInformation negativity(ρ ⊗ σ, [2, 2], 2) negativity(proj((1/sqrt(2)*(ket(0,2)⊗ket(0,2)-ket(1,2)⊗ket(1,2)))), [2, 2], 2)log_negativity(ρ ⊗ σ, [2, 2], 2)\n[*Positive partial transpose*](https://en.wikipedia.org/wiki/Peres%E2%80%93Horodecki_criterion) (the Peres–Horodecki criterion) is a necessary condition of separability of the joint state \\$\\\\rho_{AB}\\$. According PPT criterion, if \\$\\\\rho^{T_B}\\$ has non negative eigenvalues, then \\$\\\\rho_{AB}\\$ is separable.\n@repl QuantumInformation ppt(ρ ⊗ σ, [2, 2], 2) ppt(proj((1/sqrt(2)*(ket(0,2)⊗ket(0,2)-ket(1,2)⊗ket(1,2)))), [2, 2], 2) ``` –>"
},

{
    "location": "man/random/#",
    "page": "Random quantum objects",
    "title": "Random quantum objects",
    "category": "page",
    "text": ""
},

{
    "location": "man/random/#Random-quantum-objects-1",
    "page": "Random quantum objects",
    "title": "Random quantum objects",
    "category": "section",
    "text": ""
},

{
    "location": "lib/QuantumInformation/#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "lib/QuantumInformation/#Documentation-1",
    "page": "Library",
    "title": "Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "lib/QuantumInformation/#Contents-1",
    "page": "Library",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [joinpath(\"content\", f) for f in readdir(\"content\")]"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.QuantumInformation",
    "page": "Library",
    "title": "QuantumInformation.QuantumInformation",
    "category": "module",
    "text": "Main module for QuantumInformation.jl – a Julia package for numerical computation in quantum information theory.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.DynamicalMatrix",
    "page": "Library",
    "title": "QuantumInformation.DynamicalMatrix",
    "category": "type",
    "text": "T: quantum channel map.\n\nRepresentation of quantum channel by Dynamical matrix operators.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.HaarKet",
    "page": "Library",
    "title": "QuantumInformation.HaarKet",
    "category": "type",
    "text": "ϕ: vector.\n\nGenerates random ket based on ϕ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.HilbertSchmidtStates",
    "page": "Library",
    "title": "QuantumInformation.HilbertSchmidtStates",
    "category": "type",
    "text": "d: length.\n\nGenerates random ket of length d.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.IdentityChannel",
    "page": "Library",
    "title": "QuantumInformation.IdentityChannel",
    "category": "type",
    "text": "T: quantum channel map.\n\nRepresentation of identity channel.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.KrausOperators",
    "page": "Library",
    "title": "QuantumInformation.KrausOperators",
    "category": "type",
    "text": "T: quantum channel map.\n\nRepresentation of quantum channel by Kraus operators.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.Stinespring",
    "page": "Library",
    "title": "QuantumInformation.Stinespring",
    "category": "type",
    "text": "T: quantum channel map.\n\nStinespring representation of quantum channel.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.SuperOperator",
    "page": "Library",
    "title": "QuantumInformation.SuperOperator",
    "category": "type",
    "text": "T: quantum channel map.\n\nRepresentation of quantum channel by super-operator.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.SuperOperator-Union{Tuple{T}, Tuple{Function,Int64,Int64}} where T<:(AbstractArray{#s3,2} where #s3<:Number)",
    "page": "Library",
    "title": "QuantumInformation.SuperOperator",
    "category": "method",
    "text": "SuperOperator(m)\n\n\nchannel: quantum channel map.\ndim: square root of the super-operator matrix dimension.\n\nTransforms quntum channel into super-operator matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.UnitaryChannel",
    "page": "Library",
    "title": "QuantumInformation.UnitaryChannel",
    "category": "type",
    "text": "T: quantum channel map.\n\nRepresentation of unitary channel.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.applychannel-Tuple{KrausOperators{#s3} where #s3<:(AbstractArray{#s2,2} where #s2<:Number),AbstractArray{#s1,2} where #s1<:Number}",
    "page": "Library",
    "title": "QuantumInformation.applychannel",
    "category": "method",
    "text": "applychannel(Φ, ρ)\n\n\nΦ: list of vectors.\nρ: input matrix.\n\nReturn application of channel Φonρ`. Kraus representation of quantum channel Phi is a set K_i_iin I of bounded operators on mathcalH such that sum_iin I K_i^dagger K_i = mathcal1. Then Phi(rho)=sum_iin I K_i rho K_i^dagger.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.applychannel-Tuple{Stinespring{#s3} where #s3<:(AbstractArray{#s2,2} where #s2<:Number),AbstractArray{#s1,2} where #s1<:Number}",
    "page": "Library",
    "title": "QuantumInformation.applychannel",
    "category": "method",
    "text": "applychannel(Φ, ρ)\n\n\nΦ: Stinespring representation of quantum channel.\nρ: quantum state.\ndims: dimensions of registers of ρ.\n\nApplication of Stinespring representation of quantum channel into state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.applychannel-Tuple{SuperOperator{#s3} where #s3<:(AbstractArray{#s2,2} where #s2<:Number),AbstractArray{#s1,2} where #s1<:Number}",
    "page": "Library",
    "title": "QuantumInformation.applychannel",
    "category": "method",
    "text": "applychannel(Φ, ρ)\n\n\nΦ: super-operator matrix.\nρ: quantum state.\n\nApplication of super-operator matrix into state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.applychannel-Union{Tuple{T}, Tuple{DynamicalMatrix{#s5} where #s5<:(AbstractArray{#s3,2} where #s3<:Number),AbstractArray{T,2}}} where T<:Number",
    "page": "Library",
    "title": "QuantumInformation.applychannel",
    "category": "method",
    "text": "applychannel(Φ, ρ)\n\n\nΦ: dynamical matrix.\nρ: quantum state.\n\nApplication of dynamical matrix into state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.bra-Tuple{Int64,Int64}",
    "page": "Library",
    "title": "QuantumInformation.bra",
    "category": "method",
    "text": "bra(val, dim)\n\n\nval: non-zero entry - label.\ndim: length of the vector\n\nReturn Hermitian conjugate langle val = valrangle^dagger of the ket with the same label.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.bures_angle-Tuple{AbstractArray{#s11,2} where #s11<:Number,AbstractArray{#s10,2} where #s10<:Number}",
    "page": "Library",
    "title": "QuantumInformation.bures_angle",
    "category": "method",
    "text": "bures_angle(ρ, σ)\n\n\nρ: quantum state.\nσ: quantum state.\n\nReturn Bures angle between quantum states ρ and σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.bures_distance-Tuple{AbstractArray{#s11,2} where #s11<:Number,AbstractArray{#s10,2} where #s10<:Number}",
    "page": "Library",
    "title": "QuantumInformation.bures_distance",
    "category": "method",
    "text": "bures_distance(ρ, σ)\n\n\nρ: quantum state.\nσ: quantum state.\n\nReturn Bures distance between quantum states ρ and σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.concurrence-Tuple{AbstractArray{#s12,2} where #s12<:Number}",
    "page": "Library",
    "title": "QuantumInformation.concurrence",
    "category": "method",
    "text": "concurrence(ρ)\n\n\nρ: quantum state.\n\nCalculates the concurrence of a two-qubit system ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.diamond_distance-Union{Tuple{T}, Tuple{DynamicalMatrix{T},DynamicalMatrix{T}}} where T<:(AbstractArray{#s12,2} where #s12<:Number)",
    "page": "Library",
    "title": "QuantumInformation.diamond_distance",
    "category": "method",
    "text": "diamond_distance(Φ1, Φ2)\n\n\nΦ1: DynamicalMatrix\nΦ2: DynamicalMatrix\n\nReturn diamond distance between dynamical matrices Φ1 and Φ2.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.fidelity-Tuple{AbstractArray{#s11,2} where #s11<:Number,AbstractArray{#s10,2} where #s10<:Number}",
    "page": "Library",
    "title": "QuantumInformation.fidelity",
    "category": "method",
    "text": "fidelity(ρ, σ)\n\n\nρ: matrix.\nσ: matrix.\n\nReturn fidelity between matrices ρ and σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.fidelity_sqrt-Tuple{AbstractArray{#s11,2} where #s11<:Number,AbstractArray{#s10,2} where #s10<:Number}",
    "page": "Library",
    "title": "QuantumInformation.fidelity_sqrt",
    "category": "method",
    "text": "fidelity_sqrt(ρ, σ)\n\n\nρ: matrix.\nσ: matrix.\n\nReturn square root of fidelity between matrices ρ and σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.gate_fidelity-Tuple{AbstractArray{#s11,2} where #s11<:Number,AbstractArray{#s10,2} where #s10<:Number}",
    "page": "Library",
    "title": "QuantumInformation.gate_fidelity",
    "category": "method",
    "text": "gate_fidelity(U, V)\n\n\nU: quantum gate.\nV: quantum gate.\n\nReturn fidelity between gates U and V.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.grover-Tuple{Int64}",
    "page": "Library",
    "title": "QuantumInformation.grover",
    "category": "method",
    "text": "grover(dim)\n\n\nd: dimension of operator.\n\nPrepares Grover operator of dimension d.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.hadamard-Tuple{Int64}",
    "page": "Library",
    "title": "QuantumInformation.hadamard",
    "category": "method",
    "text": "hadamard(dim)\n\n\nd: dimension of operator.\n\nPrepares Hadamard operator of dimension d.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.hs_distance-Tuple{AbstractArray{#s11,2} where #s11<:Number,AbstractArray{#s10,2} where #s10<:Number}",
    "page": "Library",
    "title": "QuantumInformation.hs_distance",
    "category": "method",
    "text": "hs_distance(A, B)\n\n\nA: matrix.\nB: matrix.\n\nReturn Hilbert–Schmidt distance between matrices A and B.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.iscptp-Tuple{KrausOperators{#s11} where #s11<:(AbstractArray{#s12,2} where #s12<:Number)}",
    "page": "Library",
    "title": "QuantumInformation.iscptp",
    "category": "method",
    "text": "_\n\niscptp(Φ; atol)\n\n\nΦ: list of Kraus operators.\natol: tolerance of approximation.\n\nChecks if set of Kraus operators fulfill completness relation.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.js_divergence-Tuple{AbstractArray{#s11,2} where #s11<:Number,AbstractArray{#s10,2} where #s10<:Number}",
    "page": "Library",
    "title": "QuantumInformation.js_divergence",
    "category": "method",
    "text": "js_divergence(ρ, σ)\n\n\nρ: quantum state.\nσ: quantum state.\n\nReturn Jensen–Shannon divergence of quantum state ρ with respect to σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ket-Tuple{Int64,Int64}",
    "page": "Library",
    "title": "QuantumInformation.ket",
    "category": "method",
    "text": "ket(val, dim)\n\n\nval: non-zero entry - label.\ndim: length of the vector.\n\nReturn complex column vector valrangle of unit norm describing quantum state.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ketbra-Tuple{Int64,Int64,Int64}",
    "page": "Library",
    "title": "QuantumInformation.ketbra",
    "category": "method",
    "text": "ketbra(valk, valb, dim)\n\n\nvalk: non-zero entry - label.\nvalb: non-zero entry - label.\ndim: length of the vector\n\nReturn outer product valkranglelangle vakb of states valkrangle and valbrangle.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.kl_divergence-Tuple{AbstractArray{#s11,2} where #s11<:Number,AbstractArray{#s10,2} where #s10<:Number}",
    "page": "Library",
    "title": "QuantumInformation.kl_divergence",
    "category": "method",
    "text": "kl_divergence(ρ, σ)\n\n\nρ: quantum state.\nσ: quantum state.\n\nReturn Kullback–Leibler divergence of quantum state ρ with respect to σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.log_negativity-Tuple{AbstractArray{#s12,2} where #s12<:Number,Array{Int64,1},Int64}",
    "page": "Library",
    "title": "QuantumInformation.log_negativity",
    "category": "method",
    "text": "log_negativity(ρ, dims, sys)\n\n\nρ: quantum state.\ndims: dimensions of subsystems.\nsys: transposed subsystem.\n\nReturn log negativity of quantum state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.max_entangled-Tuple{Int64}",
    "page": "Library",
    "title": "QuantumInformation.max_entangled",
    "category": "method",
    "text": "max_entangled(d)\n\n\nd: length of the vector.\n\nReturn maximally entangled state frac1sqrtdsum_i=0^sqrtd-1iirangle of length sqrtd.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.max_mixed-Tuple{Int64}",
    "page": "Library",
    "title": "QuantumInformation.max_mixed",
    "category": "method",
    "text": "max_mixed(d)\n\n\nd: length of the vector.\n\nReturn maximally mixed state frac1dsum_i=0^d-1iranglelangle i  of length d.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.negativity-Tuple{AbstractArray{#s12,2} where #s12<:Number,Array{Int64,1},Int64}",
    "page": "Library",
    "title": "QuantumInformation.negativity",
    "category": "method",
    "text": "negativity(ρ, dims, sys)\n\n\nρ: quantum state.\ndims: dimensions of subsystems.\nsys: transposed subsystem.\n\nReturn negativity of quantum state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.norm_diamond-Union{Tuple{DynamicalMatrix{T}}, Tuple{T}} where T<:(AbstractArray{#s12,2} where #s12<:Number)",
    "page": "Library",
    "title": "QuantumInformation.norm_diamond",
    "category": "method",
    "text": "norm_diamond(Φ)\n\n\nΦ: DynamicalMatrix\n\nReturn diamond norm of dynamical matrix Φ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.norm_hs-Tuple{AbstractArray{#s12,2} where #s12<:Number}",
    "page": "Library",
    "title": "QuantumInformation.norm_hs",
    "category": "method",
    "text": "norm_hs(A)\n\n\nA: matrix.\n\nReturn Hilbert–Schmidt norm of matrix A.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.norm_trace-Tuple{AbstractArray{#s12,2} where #s12<:Number}",
    "page": "Library",
    "title": "QuantumInformation.norm_trace",
    "category": "method",
    "text": "norm_trace(A)\n\n\nA: matrix.\n\nReturn trace norm of matrix A.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.permutesystems-Union{Tuple{T}, Tuple{AbstractArray{T,2},Array{Int64,1},Array{Int64,1}}} where T<:Number",
    "page": "Library",
    "title": "QuantumInformation.permutesystems",
    "category": "method",
    "text": "permutesystems(ρ, dims, systems)\n\n\nρ: input state.\ndims: dimensions of registers of ρ.\nsystems: permuted registers.\n\nReturns state ρ with permuted registers denoted by systems.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ppt-Tuple{AbstractArray{#s12,2} where #s12<:Number,Array{Int64,1},Int64}",
    "page": "Library",
    "title": "QuantumInformation.ppt",
    "category": "method",
    "text": "ppt(ρ, dims, sys)\n\n\nρ: quantum state.\ndims: dimensions of subsystems.\nsys: transposed subsystem.\n\nReturn minimum eigenvalue of positive partial transposition of quantum state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.proj-Tuple{AbstractArray{#s2,1} where #s2<:Number}",
    "page": "Library",
    "title": "QuantumInformation.proj",
    "category": "method",
    "text": "proj(ψ)\n\n\nket: input column vector.\n\nReturn outer product ketranglelangle ket of ket.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ptrace-Tuple{AbstractArray{#s11,1} where #s11<:Number,Array{Int64,1},Int64}",
    "page": "Library",
    "title": "QuantumInformation.ptrace",
    "category": "method",
    "text": "ptrace(ψ, idims, sys)\n\n\nψ: quantum state pure state (ket).\nidims: dimensins of subsystems - only bipartite states accepted.\nsys: traced subsystem.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ptrace-Tuple{AbstractArray{#s12,2} where #s12<:Number,Array{Int64,1},Array{Int64,1}}",
    "page": "Library",
    "title": "QuantumInformation.ptrace",
    "category": "method",
    "text": "ptrace(ρ, idims, isystems)\n\n\nρ: quantum state.\nidims: dimensins of subsystems.\nisystems: traced subsystems.\n\nReturn partial trace of matrix ρ over the subsystems determined by isystems.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ptrace-Tuple{AbstractArray{#s12,2} where #s12<:Number,Array{Int64,1},Int64}",
    "page": "Library",
    "title": "QuantumInformation.ptrace",
    "category": "method",
    "text": "ptrace(ρ, idims, sys)\n\n\nρ: quantum state.\nidims: dimensins of subsystems.\nsys: traced subsystem.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ptranspose-Tuple{AbstractArray{#s11,2} where #s11<:Number,Array{Int64,1},Array{Int64,1}}",
    "page": "Library",
    "title": "QuantumInformation.ptranspose",
    "category": "method",
    "text": "ptranspose(ρ, idims, isystems)\n\n\nρ: quantum state.\nidims: dimensins of subsystems.\nisystems: transposed subsystems.\n\nReturn partial transposition of matrix ρ over the subsystems determined by isystems.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ptranspose-Tuple{AbstractArray{#s12,2} where #s12<:Number,Array{Int64,1},Int64}",
    "page": "Library",
    "title": "QuantumInformation.ptranspose",
    "category": "method",
    "text": "ptranspose(ρ, idims, sys)\n\n\nρ: quantum state.\nidims: dimensins of subsystems.\nsys: transposed subsystem.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.qft-Tuple{Int64}",
    "page": "Library",
    "title": "QuantumInformation.qft",
    "category": "method",
    "text": "qft(d)\n\n\nd: dimension of operator.\n\nPrepares gate realized a quantum Fourier transform of dimension d.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.relative_entropy-Tuple{AbstractArray{#s11,2} where #s11<:Number,AbstractArray{#s10,2} where #s10<:Number}",
    "page": "Library",
    "title": "QuantumInformation.relative_entropy",
    "category": "method",
    "text": "relative_entropy(ρ, σ)\n\n\nρ: quantum state.\nσ: quantum state.\n\nReturn quantum relative entropy of quantum state ρ with respect to σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.res-Tuple{AbstractArray{#s2,2} where #s2<:Number}",
    "page": "Library",
    "title": "QuantumInformation.res",
    "category": "method",
    "text": "res(ρ)\n\n\nρ: input matrix.\n\nReturns vec(ρ.T). Reshaping maps     matrix ρ into a vector row by row.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.reshuffle-Tuple{AbstractArray{#s11,2} where #s11<:Number}",
    "page": "Library",
    "title": "QuantumInformation.reshuffle",
    "category": "method",
    "text": "reshuffle(ρ)\n\n\nρ: reshuffled matrix.\n\nPerforms reshuffling of indices of a matrix.   Given multiindexed matrix M_(mμ)(nν) it returns   matrix M_(mn)(μν).\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.shannon_entropy-Tuple{AbstractArray{#s12,1} where #s12<:Real}",
    "page": "Library",
    "title": "QuantumInformation.shannon_entropy",
    "category": "method",
    "text": "shannon_entropy(p)\n\n\np: vector.\n\nReturn Shannon entorpy of vector p.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.shannon_entropy-Tuple{Real}",
    "page": "Library",
    "title": "QuantumInformation.shannon_entropy",
    "category": "method",
    "text": "shannon_entropy(x)\n\n\nx: real number.\n\nReturn binary Shannon entorpy given by -x  log(x) - (1 - x)  log(1 - x).\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.superfidelity-Tuple{AbstractArray{#s11,2} where #s11<:Number,AbstractArray{#s10,2} where #s10<:Number}",
    "page": "Library",
    "title": "QuantumInformation.superfidelity",
    "category": "method",
    "text": "superfidelity(ρ, σ)\n\n\nρ: quantum state.\nσ: quantum state.\n\nReturn superfidelity between quantum states ρ and σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.trace_distance-Union{Tuple{T2}, Tuple{T1}, Tuple{AbstractArray{T1,2},AbstractArray{T2,2}}} where T2<:Number where T1<:Number",
    "page": "Library",
    "title": "QuantumInformation.trace_distance",
    "category": "method",
    "text": "trace_distance(A, B)\n\n\nA: matrix.\nB: matrix.\n\nReturn trace distance between matrices A and B.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.unres-Tuple{AbstractArray{#s2,1} where #s2<:Number}",
    "page": "Library",
    "title": "QuantumInformation.unres",
    "category": "method",
    "text": "unres(ρ)\n\n\nϕ: input matrix.\n\nReturn de-reshaping of the vector into a matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.vonneumann_entropy-Tuple{LinearAlgebra.Hermitian{#s12,S} where S<:(AbstractArray{#s571,2} where #s571<:#s12) where #s12<:Number}",
    "page": "Library",
    "title": "QuantumInformation.vonneumann_entropy",
    "category": "method",
    "text": "vonneumann_entropy(ρ)\n\n\nρ: quantum state.\n\nReturn Von Neumann entropy of quantum state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.werner_state-Tuple{Int64,Float64}",
    "page": "Library",
    "title": "QuantumInformation.werner_state",
    "category": "method",
    "text": "werner_state(d, α)\n\n\nd: length of the vector.\nα: real number from [0, 1].\n\nReturns Werner state given by $ \\frac{\\alpha}{d}\\Big(\\sum{i=0}^{\\sqrt{d}-1}|ii\\rangle\\Big) \\Big(\\sum{i=0}^{\\sqrt{d}-1}\\langle ii|\\Big)+ \\frac{1-\\alpha}{d}\\sum_{i=0}^{d-1}|i\\rangle\\langle i |$.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{DynamicalMatrix{T1}},KrausOperators{T2}}} where T2<:(AbstractArray{#s3,2} where #s3<:Number) where T1<:(AbstractArray{#s5,2} where #s5<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\nΦ: list of Kraus operators.\n\nTransforms list of Kraus operators into dynamical matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{DynamicalMatrix{T1}},SuperOperator{T2}}} where T2<:(AbstractArray{#s3,2} where #s3<:Number) where T1<:(AbstractArray{#s5,2} where #s5<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\nΦ: super-operator matrix.\n\nTransforms super-operator matrix into dynamical matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{KrausOperators{T1}},DynamicalMatrix{T2}}} where T2<:(AbstractArray{#s1,2} where #s1<:Number) where T1<:(AbstractArray{#s2,2} where #s2<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\nΦ: dynamical matrix.\n\nTransforms dynamical matrix into list of Kraus operators.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{KrausOperators{T1}},SuperOperator{T2}}} where T2<:(AbstractArray{#s3,2} where #s3<:Number) where T1<:(AbstractArray{#s5,2} where #s5<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\nΦ: super-operator matrix.\n\nTransforms super-operator matrix into list of Kraus operators.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{Stinespring{T1}},DynamicalMatrix{T2}}} where T2<:(AbstractArray{#s3,2} where #s3<:Number) where T1<:(AbstractArray{#s5,2} where #s5<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\nΦ: dynamical matrix.\n\nTransforms dynamical matrix into Stinespring representation of quantum channel.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{Stinespring{T1}},KrausOperators{T2}}} where T2<:(AbstractArray{#s1,2} where #s1<:Number) where T1<:(AbstractArray{#s2,2} where #s2<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\nΦ: list of Kraus operators.\n\nTransforms list of Kraus operators into Stinespring representation of quantum channel.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{Stinespring{T1}},SuperOperator{T2}}} where T2<:(AbstractArray{#s3,2} where #s3<:Number) where T1<:(AbstractArray{#s5,2} where #s5<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\nΦ: super-operator matrix.\n\nTransforms super-operator matrix into Stinespring representation of quantum channel.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{SuperOperator{T1}},DynamicalMatrix{T2}}} where T2<:(AbstractArray{#s3,2} where #s3<:Number) where T1<:(AbstractArray{#s5,2} where #s5<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\nΦ: dynamical matrix.\n\nTransforms dynamical matrix into super-operator matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{SuperOperator{T1}},KrausOperators{T2}}} where T2<:(AbstractArray{#s3,2} where #s3<:Number) where T1<:(AbstractArray{#s5,2} where #s5<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\nΦ: list of Kraus operators.\n\nTransforms list of Kraus operators into super-operator matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.purity-Tuple{AbstractArray{#s12,2} where #s12<:Number}",
    "page": "Library",
    "title": "QuantumInformation.purity",
    "category": "method",
    "text": "purity(ρ)\n\n\nρ: matrix.\n\nReturn the purity of ρ ∈ [1/d, 1]\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Index-1",
    "page": "Library",
    "title": "Index",
    "category": "section",
    "text": "A list of all documentation sorted by module.Modules = [QuantumInformation]"
},

]}

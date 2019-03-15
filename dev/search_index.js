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
    "text": "A Julia package for numerical computation in quantum information theory.Numerical investigations are prevalent in quantum information theory. Numerical experiments can be used to find counter examples for theorems, to test hypotheses or to gain insight about quantum objects and operations.Our goal while designing QuantumInformation.jl library was to follow principles presented in book \"Geometry of Quantum States\" [1]. We work with column vectors reprinting kets and row vectors representing bras. We fix our basis to the computational one. Density matrices and quantum channels are represented as two dimensional arrays in the same fixed basis. This approach allows us to obtain low level complexity of our code, high flexibility and good computational efficiency. The design choices where highly motivated by the properties of the language in which the our library was implemented, namely Julia [2]."
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
    "page": "Linear algebra in Julia",
    "title": "Linear algebra in Julia",
    "category": "page",
    "text": ""
},

{
    "location": "man/vectors/#Linear-algebra-in-Julia-1",
    "page": "Linear algebra in Julia",
    "title": "Linear algebra in Julia",
    "category": "section",
    "text": "A basic construction of vector in Julia creates a full one-index array containing elements of a number type as presented below.julia> x = [0.0, 1.0im]\n2-element Array{Complex{Float64},1}:\n 0.0+0.0im\n 0.0+1.0imA transposition of a column vector return an object of type LinearAlgebra.Transpose as shown belowjulia> xt = transpose(x)\n1×2 LinearAlgebra.Transpose{Complex{Float64},Array{Complex{Float64},1}}:\n 0.0+0.0im  0.0+1.0imWhile a~Hermitian conjugate of the same vector returns a LinearAlgebra.Adjoint parametrized by the type Array:julia> xc = [0.0, 1.0im]\'\n1×2 LinearAlgebra.Adjoint{Complex{Float64},Array{Complex{Float64},1}}:\n 0.0-0.0im  0.0-1.0imValues of variables xt and xc are views of the value of variable x. The column and row vectors behave like bras and kets, for example xc*x denotes the inner product of bra xc and ket x, while x*xc denotes its outer product resulting in a two-index array.The linear algebra library in Julia provides standard operations on matrices and vectors that are designed to take in to the account the types of objects."
},

{
    "location": "man/states/#",
    "page": "States and channels",
    "title": "States and channels",
    "category": "page",
    "text": "using QuantumInformation"
},

{
    "location": "man/states/#States-and-channels-1",
    "page": "States and channels",
    "title": "States and channels",
    "category": "section",
    "text": "In this and the following sections we will denote complex Euclidean spaces mathbbC^d with mathcalX, mathcalY, mathcalZ etc. When needed the dimension of a space mathcalX will be denoted mathrmdim(mathcalX). The set of matrices transforming vectors from mathcalX to mathcalY will be denoted mathrmL(mathcalX mathcalY). For simplicity we will write mathrmL(mathcalX) equiv mathrmL(mathcalX mathcalX)."
},

{
    "location": "man/states/#States-1",
    "page": "States and channels",
    "title": "States",
    "category": "section",
    "text": "By psirangleinmathcalX we denote a normed column vector. Notice that any psirangle can be expressed as psirangle=sum_i=1^n alpha_i irangle, where sum_i=1^n alpha_i^2=1 and the set irangle_i=1^n is the computational basis.ket(1,2)\n(1/sqrt(2)) * (ket(1,2) + ket(2,2))According to common academic convention, we count the indices of states starting from one. Following the standard Dirac notation the symbol langlepsi denotes the row vector dual to psirangle. Therefore psirangle=langlepsi^dagger, where the symbol ^dagger denotes the Hermitian conjugation.bra(2,3)  The inner product of phirangle psirangle in mathcalX is denoted by langlepsiphirangle and the norm is defined as phirangle=sqrtlanglephiphirangle.ψ=(1/sqrt(2)) * (ket(1,2) + ket(2,2))\nϕ=(1/2) * ket(1,2) + (sqrt(3)/2) * ket(2,2)\n\nϕ\'*ψ\n\nsqrt(ϕ\'*ϕ)The form psiranglelanglephi denotes outer product of psirangle and langlephi from mathrmL(mathcalX).ketbra(2,3,4)Specifically, psiranglelanglepsi is a rank-one projection operator called as pure state. Generally, any quantum state rho can be expressed as rho=sum_i=0^n q_i iranglelangle i, where sum_i=0^n q_i=1. Notice that rho is a trace-one positive semi-definite linear operator i.e.: rho=rho^dagger, rhogeq 0 and mathrmtrrho=1.proj(ψ)For convenience, the QuantumInformation.jl library provides the implementations of maximally mixed, maximally entangled and Werner states.max_entangled(4)\nmax_mixed(4)\nwerner_state(4, 0.4)"
},

{
    "location": "man/states/#Non-standard-matrix-transformations-1",
    "page": "States and channels",
    "title": "Non-standard matrix transformations",
    "category": "section",
    "text": "We will now introduce reshaping operators, which map matrices to vectors and vice versa. We start with the mapping mathrmresmathrmL(mathcalXY)tomathcalYotimesmathcalX, which transforms the matrix rho into a vector row by row. More precisely, for dyadic operators psiranglelanglephi, where psirangle in mathcalY, phirangle in mathcalX the operation mathrmres is defined as mathrmres(psiranglelanglephi)=psirangleoverlinephirangle and can be uniquely extend to the whole space mathrmL(mathcalXY) by linearity.res(ketbra(1,2,2))The inverse operation to mathrmres is mathrmunresmathcalYotimesmathcalXto mathrmL(mathcalXY)  which transforms the vector into a matrix It is defined as the unique linear mapping satisfying rho=mathrmunres(mathrmres(rho)).unres(res(ketbra(1,2,2)))Let us recall that trace is a mapping mathrmTrmathrmL(mathcalX)to mathbbC given by mathrmTrrhomapstosum_i=1^mathrmdim(mathcalX)langle e_irhoe_irangle, where e_irangle  is an orthonormal basis of mathcalX. According to this, partial trace is a mapping mathrmTr_mathcalX mathrmL(mathcalXotimesmathcalY) to mathrmL(mathcalY) such that mathrmTr_mathcalX rho_Aotimes rho_B mapsto rho_B mathrmTr(rho_A), where rho_Ain mathrmL(mathcalX), rho_Bin mathrmL(mathcalY). As this is a linear map, it may be uniquely extended to the case of operators which are not in a tensor product form.ρ = [0.25 0.25im; -0.25im 0.75]\nσ = [0.4 0.1im; -0.1im 0.6]\n\nptrace(ρ ⊗ σ, [2, 2], [2])Matrix transposition is a mapping ^TmathrmL(mathcalXY) to mathrmL(mathcalYX) such that left(rho^T right)_ij = rho_ji, where rho_ij is a i-th row, j-th column element of matrix rho. Following this, we may introduce \\emph{partial transposition} ^Gamma_B mathrmL(mathcalX_A otimes mathcalX_B mathcalY_A otimes mathcalY_B) to mathrmL(mathcalX_A otimes mathcalY_B mathcalY_A otimes mathcalX_B), which for a product state rho_Aotimesrho_B is given by ^Gamma_B rho_Aotimesrho_Bmapstorho_Aotimesrho_B^T. The definition of partial transposition can be uniquely extended for all operators from linearity.ptranspose(ρ ⊗ σ, [2, 2], [1])For given multiindexed matrix rho_(mmu)(nnu)=langle m murhon nurangle, the reshuffle operation is defined as rho^R_(mmu)(nnu)=rho_(mn)(munu).reshuffle(ρ ⊗ σ)"
},

{
    "location": "man/states/#Channels-1",
    "page": "States and channels",
    "title": "Channels",
    "category": "section",
    "text": "Physical transformations of quantum states into quantum states are called quantum channels i.e. linear Completely Positive Trace Preserving (CP-TP) transformations. Probabilistic transformations of quantum states are called quantum operations and mathematically they are defined as linear Completely Positive Trace Non-increasing (CP-TNI) maps. For the sake of simplicity we will refer to both CP-TP and CP-TNI maps as quantum channels when it will not cause confusion.There exists various representations of quantum channels such as:Kraus operators,\nnatural representation, also called superoperator representation,\nStinespring representation,\nChoi-Jamiołkowski matrices, sometimes called dynamical matrices.The product of superoperators Phi_1in mathrmT(mathcalX_1mathcalY_1), Phi_2in mathrmT(mathcalX_2mathcalY_2) is a mapping Phi_1otimesPhi_2in T(mathcalX_1otimesmathcalX_2mathcalY_1otimesmathcalY_2) that satisfies (Phi_1otimesPhi_2)(rho_1otimesrho_2)=Phi_1(rho_1)otimesPhi_2(rho_2). For the operators that are not in a tensor product form this notion can be uniquely extended from linearity.According to Kraus\' theorem, any completely positive trace-preserving (CP-TP) map Phi can always be written as Phi(rho)=sum_i=1^r K_i rho K_i^dagger for some set of operators K_i_i=1^r satisfying sum_i=1^r K_i^dagger K_i = mathbbI_mathcalX, where r is the rank of superoperator Phi.Another way to represent the quantum channel is based on Choi-Jamiołkowski isomorphism. Consider mapping JmathrmT(mathcalXY)to mathrmL(mathcalYotimesmathcalX) such that J(Phi)=(PhiotimesmathbbI_mathrmL(mathcalX)) (mathrmres(mathbbI_mathcalX) mathrmres(mathbbI_mathcalX)^dagger). Equivalently J(Phi)=sum_ij=1^mathrmdim(mathcalX)Phi(iranglelangle j)otimesiranglelangle j. The action of a superoperator in the Choi representation is given by Phi(rho)=mathrmTr_mathcalX(J(Phi)(mathbbI_mathcalYotimesrho^T)).The natural representation of a quantum channel mathrmT(mathcalX mathcalY) is a mapping mathrmres(rho) mapsto mathrmres(Phi(rho)). It is represented by a matrix K(Phi) in mathrmL(mathcalX otimes mathcalX mathcalY otimes mathcalY) for which the following holds \\begin{equation} K(\\Phi) \\mathrm{res}(\\rho) = \\mathrm{res}(\\Phi(\\rho)), \\end{equation} for all rho in mathrmL(mathcalX).Let mathcalX mathcalY and mathcalZ be a complex Euclidean spaces. The action of the Stinespring representation of a quantum channel Phiin mathrmT(mathcalXmathcalY) on a state rhoin mathrmL(mathcalX) is given by \\begin{equation} \\Phi(\\rho)=\\mathrm{Tr}_\\mathcal{Z}(A\\rho A^\\dagger), \\end{equation} where AinmathrmL(mathcalXmathcalYotimesmathcalZ).We now briefly describe the relationships among channel representations [1]. Let Phiin mathrmT(mathcalX mathcalY) be a quantum channel which can be written in the Kraus representation as Phi(rho)=sum_i=1^r K_i rho K_i^dagger, where K_i_i=1^r are Kraus operators satisfying sum_i=1^r K_i^dagger K_i = mathbbI_mathcalX. According to this assumption, Phi can be represented inChoi representation as J(Phi)=sum_i=1^r mathrmres(K_i)mathrmres(K_i^dagger),\nnatural representation as K(Phi)=sum_i=1^r K_iotimes K_i^*,\nStinespring representation as Phi(rho)=mathrmTr_mathcalZ(Arho A^dagger),where A=sum_i=1^r K_iotimes e_irangle and mathcalZ=mathbbC^r.In QuantumInformation.jl states and channels are always represented in the computational basis therefore channels are stored in the memory as either vectors of matrices in case of Kraus operators or matrices in other cases. In QuantumInformation.jl quantum channels are represented by a set of types deriving from an abstract type AbstractQuantumOperation{T} where type parameter T should inherit from AbstractMatrix{<:Number}. Every type inheriting from AbstractQuantumOperation{T} should contain fields idim and odim representing the dimension of input and output space of the quantum channel.Two special types of channels are implemented: UnitaryChannel and IdentityChannel that can transform ket vectors into ket vectors."
},

{
    "location": "man/states/#Constructors-1",
    "page": "States and channels",
    "title": "Constructors",
    "category": "section",
    "text": "Channel objects can be constructed from matrices that represent them, as shown in the following listingγ=0.4\nK0 = Matrix([1 0; 0 sqrt(1-γ)])\nK1 = Matrix([0 sqrt(γ); 0 0])\n\nΦ = KrausOperators([K0,K1])\n\niscptp(Φ)There are no checks whether a matrix represents a valid CP-TP or CP-TNI map, because this kind of verification is costly and requires potentially expensive numerical computation. Function such as iscptp(), and iscptni() are provided to test properties of supposed quantum channel or quantum operation."
},

{
    "location": "man/states/#Conversion-1",
    "page": "States and channels",
    "title": "Conversion",
    "category": "section",
    "text": "Conversions between all quantum channel types, i.e. these that derive from AbstractQuantumOperation{T} are implemented. The users are not limited by any single channel representation and can transform between representations they find the most efficient or suitable for their purpose.Ψ1 = convert(SuperOperator{Matrix{ComplexF64}}, Φ)\nΨ2 = convert(DynamicalMatrix{Matrix{Float64}}, Φ)\nΨ3 = convert(Stinespring{Matrix{Float64}}, Φ)"
},

{
    "location": "man/states/#Application-1",
    "page": "States and channels",
    "title": "Application",
    "category": "section",
    "text": "Channels can act on pure and mixed states represented by vectors and matrices respectively. Channels are callable and therefore mimic application of a function on a quantum state.ρ1=ψ * ψ\'\nΦ(ρ1)\nΨ1(ρ1)\nΦ(ψ)"
},

{
    "location": "man/states/#Composition-1",
    "page": "States and channels",
    "title": "Composition",
    "category": "section",
    "text": "Channels can be composed in parallel or in sequence. Composition in parallel is done using kron() function or the overloaded otimes operator. Composition in sequence can be done in two ways either by using Julia built-in function composition operator (fcirc g)(cdot)=f(g)(cdot) or by using multiplication of objects inheriting from AbstractQuantumOperation{T} abstract type.ρ2=ϕ * ϕ\'\n\n(Φ⊗Φ)(ρ1⊗ρ2)\n(Ψ1∘Ψ2)(ρ1)"
},

{
    "location": "man/states/#refs_sc-1",
    "page": "States and channels",
    "title": "References",
    "category": "section",
    "text": "[1] J. Watrous, The Theory of Quantum Information, Cambridge University Press (2018)."
},

{
    "location": "man/functionals/#",
    "page": "Functionals",
    "title": "Functionals",
    "category": "page",
    "text": "using QuantumInformation\nusing LinearAlgebra"
},

{
    "location": "man/functionals/#Functionals-1",
    "page": "Functionals",
    "title": "Functionals",
    "category": "section",
    "text": ""
},

{
    "location": "man/functionals/#Trace-norm-and-distance-1",
    "page": "Functionals",
    "title": "Trace norm and distance",
    "category": "section",
    "text": "Let rho sigma in mathrmL(mathcalX). The trace norm is defined as rho_1 = mathrmTr sqrtrhorho^dagger and the trace distance is defined as D_1(rhosigma)=frac12rho-sigma_1.ψ=(1/sqrt(2)) * (ket(1,2) + ket(2,2))\nϕ=(1/2) * ket(1,2) + (sqrt(3)/2) * ket(2,2)\n\nρ=proj(ψ)\nσ=proj(ϕ)\n\nnorm_trace(ρ)\ntrace_distance(ρ, σ)"
},

{
    "location": "man/functionals/#Hilbert–Schmidt-norm-and-distance-1",
    "page": "Functionals",
    "title": "Hilbert–Schmidt norm and distance",
    "category": "section",
    "text": "The Hilbert–Schmidt norm norm and distance defined by rho_HS=sqrtmathrmTrrho^dagger rho and D_HS(rhosigma)=frac12rho-sigma_HS, respectively, can be used as followsnorm_hs(ρ)\nhs_distance(ρ, σ)"
},

{
    "location": "man/functionals/#Fidelity-and-superfidelity-1",
    "page": "Functionals",
    "title": "Fidelity and superfidelity",
    "category": "section",
    "text": "Fidelity is a measure of distance of quantum states. It is an example of a distance measure which is not a metric on the space of quantum states. The fidelity of two quantum states rho sigma in mathrmL(mathcalX) is given by F(rhosigma)=sqrtrhosqrtsigma_1fidelity_sqrt(ρ, σ)\nfidelity(ρ, σ)\nfidelity(ψ, σ)\nfidelity(ρ, ϕ)\nfidelity(ψ, ϕ)"
},

{
    "location": "man/functionals/#Superfidelity-1",
    "page": "Functionals",
    "title": "Superfidelity",
    "category": "section",
    "text": "Superfidelity is an upper bound on the fidelity of two quantum states It is defined by G(rho sigma) = mathrmTrrho sigma + sqrt1 - mathrmTrrho^2 sqrt1-mathrmTr sigma^2.superfidelity(ρ, σ)"
},

{
    "location": "man/functionals/#Diamond-norm-1",
    "page": "Functionals",
    "title": "Diamond norm",
    "category": "section",
    "text": "In order to introduce the \\emph{diamond norm}, we first introduce the notion of the induced trace norm. Given Phi in mathrmT(mathcalX mathcalY) we define its induced trace norm as  Phi _1 = mathrmmax left  Phi(X) _1 X in L(mathcalX)  X _1 leq 1 right. The diamond norm of Phi is defined as  Phi _diamond =  Phi otimes mathbbI_mathrmL(mathcalY) _1 One important property of the diamond norm is that for Hermiticity-preserving Phi in mathrmT(mathcalX mathcalY) we obtain  Phi _diamond = max left left (Phi otimes mathbbI_mathrmL(mathcalY)) 	left(psiranglelanglepsi right )right_1 psirangle in 	mathcalX otimes 	mathcalY 	langlepsipsirangle=1 right.γ = 0.3\nK0 = Matrix([1 0; 0 sqrt(1-γ)])\nK1 = Matrix([0 sqrt(γ); 0 0])\n\nΦ = convert(DynamicalMatrix{Array{Float64,2}}, KrausOperators([K0,K1]))\n\nL0 = Matrix([1 0; 0 sqrt(1-γ)])\nL1 = Matrix([0 0; 0 sqrt(γ)])\n\nΨ = convert(DynamicalMatrix{Array{Float64,2}}, KrausOperators([K0,K1]))\n\nnorm_diamond(Φ)\ndiamond_distance(Φ, Ψ)"
},

{
    "location": "man/functionals/#Shannon-entropy-and-von-Neumann-entropy-1",
    "page": "Functionals",
    "title": "Shannon entropy and von Neumann entropy",
    "category": "section",
    "text": "Shannon entropy is defined for a probability vector p as H(mathrmp)=-sum_i=1^n p_ilog_2 p_i. We also provide an implementation for the point Shannon entropy. It is defined as h(a) = -a log a - (1-a)log(1-a).p = vec([0.3 0.2 05])\n\nshannon_entropy(p)\nshannon_entropy(0.5)For a quantum system described by a state $\\rho$, the von Neumann entropy is S(rho)=-mathrmTr rho log rho. Let lambda_i,  0leq i  n be the eigenvalues of rho, then S(rho) can be written as S(rho)=-sum_i=1^n lambda_i log lambda_i.ρ = [0.25 0.25im; -0.25im 0.75]\nσ = [0.4 0.1im; -0.1im 0.6]\n\nvonneumann_entropy(0.4 * ρ + 0.6 * σ)"
},

{
    "location": "man/functionals/#Distinguishability-between-two-quantum-states-1",
    "page": "Functionals",
    "title": "Distinguishability between two quantum states",
    "category": "section",
    "text": "One of the measure of distinguishability between two quantum states is a qauntum relative entropy, called also Kullback–Leibler divergence, defined as S(rhosigma)=-mathrmTrrhologsigma + mathrmTrrhologrhorelative_entropy(ρ, σ)\nkl_divergence(ρ, σ)Another type of measure of distinguishability between two quantum state is quantum Jensen–Shannon divergence given by QJS(rhosigma)=Sleft(frac12rho+frac12sigmaright)- left(frac12S(rho)+frac12S(sigma)right).js_divergence(ρ, σ)The Bures distance defines an infinitesimal distance between quantum states, and it is defined as D_B=sqrt2(1-sqrtF(rhosigma)). The value related with Bures distance is the Bures angle D_A(rhosigma)=arccos(sqrtF(rhosigma))bures_distance(ρ, σ)\nbures_angle(ρ, σ)"
},

{
    "location": "man/functionals/#Quantum-entanglement-1",
    "page": "Functionals",
    "title": "Quantum entanglement",
    "category": "section",
    "text": "One of the entanglement measure is negativity defined as mathrmN(rho)=fracrho^T_A_1-12.negativity(ρ ⊗ σ, [2, 2], 2)\nnegativity(proj((1/sqrt(2)*(ket(1,2) ⊗ ket(1,2)-ket(2,2) ⊗ ket(2,2)))), [2, 2], 2)\n\nlog_negativity(ρ ⊗ σ, [2, 2], 2)Positive partial transpose (the Peres–Horodecki criterion) is a necessary condition of separability of the joint state rho_AB. According PPT criterion, if rho^T_B has non negative eigenvalues, then rho_AB is separable.ppt(ρ ⊗ σ, [2, 2], 2)\nppt(proj((1/sqrt(2)*(ket(1,2) ⊗ ket(1,2)-ket(2,2) ⊗ ket(2,2)))), [2, 2], 2)Another way to quantification of quantum entanglement is Concurrence. Concurrence of quantum state rho is a strong separability criterion. For two-qubit systems it is defined as C(rho)=max(0lambda_1-lambda_2-lambda_3-lambda_4), where lambda_i are decreasing eigenvalues of sqrtsqrtrhotilderhosqrtrho with tilderho=(sigma_yotimessigma_y)rho^*(sigma_yotimessigma_y). If C(rho)=0, then rho is separable.ρ = [0.25 0.1im; -0.1im 0.75]\nσ = [0.4 0.1im; -0.1im 0.6]\nconcurrence(ρ ⊗ σ)\nconcurrence(proj(max_entangled(4)))"
},

{
    "location": "man/measurement/#",
    "page": "Measurement",
    "title": "Measurement",
    "category": "page",
    "text": "using QuantumInformation"
},

{
    "location": "man/measurement/#Measurement-1",
    "page": "Measurement",
    "title": "Measurement",
    "category": "section",
    "text": "Measurement is modeled in two ways:as Positive Operator Valued Measures (POVMs),\nmeasurements with post-selection.In both cases a measurement is treated as a special case of a quantum channel (operation)."
},

{
    "location": "man/measurement/#Positive-Operator-Valued-Measure-measurement-1",
    "page": "Measurement",
    "title": "Positive Operator Valued Measure measurement",
    "category": "section",
    "text": "A POVM measurement is defined as follows. Let muGammatomathrmP(mathcalX) be a mapping from a finite alphabet of measurement outcomes to the set of linear positive operators. If sum_xiinGamma mu(xi)=mathbbI_mathcalX then mu is a POVM measurement. The set of positive semi-definite linear operators is defined as mathrmP(mathcalX)=Xin mathrmL(mathcalX) langlepsiXpsiranglegeq 0 text for all  psirangleinmathcalX. POVM measurement models the situation where a quantum object is destroyed during the measurement process and quantum state after the measurement does not exists.We model POVM measurement as a channel thetamathrmL(mathcalX)to mathrmL(mathcalY), where mathcalY=mathrmspanxirangle_xiinGamma such that theta(rho) = sum_xiinGamma mathrmTr(rho mu(xi))xiranglelanglexi. This channel transforms the measured quantum state into a classical state (diagonal matrix) containing probabilities of measuring given outcomes. Note that in QuantumInformation.jl Gamma=12ldotsGamma and POVM measurements are represented by the type POVMMeasurement{T} <: AbstractQuantumOperation{T} where T<:AbstractMatrix{<:Number}. Predicate function ispovm() verifies whether a list of matrices is a proper POVM.ρ=proj(1.0/sqrt(2)*(ket(1,3)+ket(3,3)))\nE0 = proj(ket(1,3))\nE1 = proj(ket(2,3))+proj(ket(3,3))\n\nM = POVMMeasurement([E0,E1])\n\nispovm(M)\nM(ρ)"
},

{
    "location": "man/measurement/#Measurement-with-post-selection-1",
    "page": "Measurement",
    "title": "Measurement with post-selection",
    "category": "section",
    "text": "When a quantum system after being measured is not destroyed one can be interested in its state after the measurement. This state depends on the measurement outcome. In this case the measurement process is defined in the following way.Let muGammato mathrmL(mathcalX mathcalY) be a mapping from a finite set of measurement outcomes to set of linear operators called effects. If sum_xiinGamma mu(xi)^dagger mu(xi)=mathbbI_mathcalX then mu is a quantum measurement. Given outcome xi was obtained, the state before the measurement, rho, is transformed into sub-normalized quantum state rho_xi=mu(xi)rhomu(xi)^dagger. The outcome xi will be obtained with probability mathrmTr(rho_xi).PM = PostSelectionMeasurement(E1)\niseffect(PM)\nPM(ρ)In QuantumInformation this kind of measurement is modeled as CP-TNI map with single Kraus operator mu(xi) and represented as PostSelectionMeasurement{T} <: AbstractQuantumOperation{T} where T<:AbstractMatrix{<:Number}. Measurement types can be composed and converted to Kraus operators, superoperators, Stinespring representation operators, and dynamical matrices.α = 0.3\nK0 = ComplexF64[0 0 sqrt(α); 0 1 0; 0 0 0]\nK1 = ComplexF64[1 0 0; 0 0 0; 0 0 sqrt(1 - α)]\nΦ = KrausOperators([K0,K1])\n\nρ=proj(1.0/sqrt(2)*(ket(1,3)+ket(3,3)))\n\n(PM∘Φ)(ρ)"
},

{
    "location": "man/random/#",
    "page": "Random quantum objects",
    "title": "Random quantum objects",
    "category": "page",
    "text": "using QuantumInformation\nusing LinearAlgebra"
},

{
    "location": "man/random/#Random-quantum-objects-1",
    "page": "Random quantum objects",
    "title": "Random quantum objects",
    "category": "section",
    "text": "In this section we present the implementation of the sub-package RandomMatrices. The justification for including these functionalities in our package is twofold. First, the application of random matrix theory (RMT) in quantum information is a blooming field of research with a plethora of interesting results [1], [2], [3], [4], [5], [6], [7], [8], [9], [10]. Hence, it is useful to have readily available implementations of known algorithms of generating random matrices. Secondly, when performing numerical investigations, we often need \'\'generic\'\' inputs. Generating random matrices with a known distribution is one of the ways to obtain such generic inputs."
},

{
    "location": "man/random/#Ginibre-matrices-1",
    "page": "Random quantum objects",
    "title": "Ginibre matrices",
    "category": "section",
    "text": "In this section we introduce the Ginibre random matrices ensemble [11]. This ensemble is at the core of a vast majority of algorithms for generating random matrices presented in later subsections. Let (G_ij)_1 leq i leq m 1 leq j leq n be a mtimes n table of independent identically distributed (i.i.d.) random variable on mathbbK. The field mathbbK can be either of mathbbR, mathbbC or mathbbQ. With each of the fields we associate a Dyson index beta equal to 1, 2, or 4 respectively. Let G_ij be i.i.d random variables with the real and imaginary parts sampled independently from the distribution mathcalN(0 frac1beta). Hence, G in mathrmL(mathcalX mathcalY), where matrix G is \\begin{equation} P(G) \\propto \\exp(-\\mathrm{Tr} G G^\\dagger). \\end{equation} This law is unitarily invariant, meaning that for any unitary matrices U and V, G and UGV are equally distributed. It can be shown that for beta=2 the eigenvalues of G are uniformly distributed over the unit disk on the complex plane.In our library the ensemble Ginibre matrices is implemented in the GinibreEnsemble{β} parametric type. The parameter determines the Dyson index. The following constructors are provided GinibreEnsemble{β}(m::Int, n::Int), GinibreEnsemble{β}(m::Int), GinibreEnsemble(m::Int, n::Int), GinibreEnsemble(m::Int). The parameters n and m determine the dimensions of output and input spaces. The versions with one argument assume m=n. When the Dyson index is omitted it assumed that beta=2. Sampling from these distributions can be performed as followsg = GinibreEnsemble{2}(2,3)\n\nrand(g)The function rand has specialized methods for each possible value of the Dyson index beta."
},

{
    "location": "man/random/#Wishart-matrices-1",
    "page": "Random quantum objects",
    "title": "Wishart matrices",
    "category": "section",
    "text": "Wishart matrices form an ensemble of random positive semidefinite matrices. They are parametrized by two factors. First is the Dyson index beta which is equal to one for real matrices, two for complex matrices and four for symplectic matrices. The second parameter, K, is responsible for the rank of the matrices. They are sampled as followsChoose beta and K.\nSample a Ginibre matrix Gin mathrmL(mathrmX mathrmY) with the Dyson index betaand mathrmdim(mathrmX) = d and mathrmdim(mathrmY)=Kd.Return GG^dagger.In QuantumInformation.jl this is implemented using the type WishartEnsemble{β, K}. We also provide additional constructors for convenience WishartEnsemble{β}(d::Int) where β = WishartEnsemble{β, 1}(d), WishartEnsemble(d::Int) = WishartEnsemble{2}(d).These can be used in the following wayw = WishartEnsemble{1,0.2}(5)\n\nz = rand(w)\n\neigvals(z)\n\nw = WishartEnsemble(3)\n\nz = rand(w)\n\neigvals(z)"
},

{
    "location": "man/random/#Circular-ensembles-1",
    "page": "Random quantum objects",
    "title": "Circular ensembles",
    "category": "section",
    "text": "Circular ensembles are measures on the space of unitary matrices. There are three main circular ensembles. Each of this ensembles has an associated Dyson index beta [12]Circular orthogonal ensemble (COE), beta=1.\nCircular unitary ensemble (CUE), beta=2.\nCircular symplectic ensemble (CSE), beta=4.They can be characterized as follows. The CUE is simply the Haar measure on the unitary group. Now, if U is an element of CUE then U^TU is an element of COE and U_R U is an element CSE. As can be seen the sampling of Haar unitaries is at the core of sampling these ensembles. Hence, we will focus on them in the remainder of this section.There are several possible approaches to generating random unitary matrices according to the Haar measure. One way is to consider known parametrizations of unitary matrices, such as the Euler [13] or Jarlskog [14] ones. Sampling these parameters from appropriate distributions yields a Haar random unitary. The downside is the long computation time, especially for large matrices, as this involves a lot of matrix multiplications. We will not go into this further, we refer the interested reader to the papers on these parametrizations.Another approach is to consider a Ginibre matrix G in mathrmL(mathrmX) and its polar decomposition G=U P, where U in mathrmL(mathrmX) is unitary and P is a positive matrix. The matrix P is unique and given by sqrtG^dagger G. Hence, assuming P is invertible, we could recover U as \\begin{equation} U = G (G^\\dagger G) ^{-\\frac{1}{2}}. \\end{equation} As this involves the inverse square root of a matrix, this approach can be potentially numerically unstable.The optimal approach is to utilize the QR decomposition of G, G=QR, where Q in mathrmL(mathrmX) is unitary and R in mathrmL(mathrmX) is upper triangular. This procedure is unique if G is invertible and we require the diagonal elements of R to be positive. As typical implementations of the QR algorithm do not consider this restriction, we must enforce it ourselves. The algorithm is as followsGenerate a Ginibre matrix G in mathrmL(mathrmX), mathrmdim(mathrmX) = d withDyson index beta=2.Perform the QR decomposition obtaining Q and R.\nMultiply the i\\textsuperscript{th} column of Q by r_iir_ii.This gives us a Haar distributed random unitary. This procedure can be generalized in order to obtain a random isometry. The only required changed is the dimension of G. We simply start with G in mathrmL(mathrmX mathrmY), where dim(mathrmX)geq dim(mathrmY).Furthermore, we may introduce two additional circular ensembles corresponding to the Haar measure on the orthogonal and symplectic groups. These are the circular real ensemble (CRE) and circular quaternion ensemble (CQE). Their sampling is similar to sampling from CUE. The only difference is the initial Dyson index of the Ginibre matrix. This is set to beta=1 for CRE and beta=4 for CQE.In QuantumInformation.jl these distributions can be sampled asc = CircularEnsemble{2}(3)\n\nu = rand(c)\n\nu*u\'Sampling from the Haar measure on the orthogonal group can be achieved asc = CircularRealEnsemble(3)\n\no = rand(c)\n\no*o\'For convenience we provide the following type aliases const COE = CircularEnsemble{1}, const CUE = CircularEnsemble{2}, const CSE = CircularEnsemble{4}."
},

{
    "location": "man/random/#Random-quantum-states-1",
    "page": "Random quantum objects",
    "title": "Random quantum states",
    "category": "section",
    "text": ""
},

{
    "location": "man/random/#Pure-states-1",
    "page": "Random quantum objects",
    "title": "Pure states",
    "category": "section",
    "text": "Pure states are elements of the unit sphere in mathrmX. Thus it is straightforward to generate them randomly. We start with a vector of dim(mathrmX) independent complex numbers sampled from the standard normal distribution. What remains is to normalize the length of this vector to unity.This is implemented using the HaarKet{β} type. The value beta=1 corresponds to the Haar measure on the unit sphere in mathbbR^d, while beta=2 corresponds to the Haar measure on the unit sphere in mathbbC^d. The usage is as followsh = HaarKet{2}(3)\n\nψ = rand(h)\n\nnorm(ψ)For convenience we provide the following constructor HaarKet(d::Int) = HaarKet{2}(d) as the majority of uses cases require sampling complex states."
},

{
    "location": "man/random/#Mixed-states-1",
    "page": "Random quantum objects",
    "title": "Mixed states",
    "category": "section",
    "text": "Random mixed states can be generated in one of two equivalent ways. The first one comes from the partial trace of random pure states. Suppose we have a pure state psirangle in mathrmX otimes mathrmY. Then we can obtain a random mixed as \\begin{equation} \\rho = \\mathrm{Tr}_\\mathrm{Y} |\\psi\\rangle\\langle\\psi|. \\end{equation} Note that in the case dim(mathrmX)=dim(mathrmY) we recover the (flat) Hilbert-Schmidt distribution on the set of quantum states.An alternative approach is to start with a Ginibre matrix G in mathrmL(mathrmX mathrmY). We obtain a random quantum state rho as \\begin{equation} \\rho = GG^\\dagger/\\mathrm{Tr}(GG^\\dagger). \\end{equation} It can be easily verified that this approach is equivalent to the one utilizing random pure states. First, note that in both cases we start with dim(mathrmX) dim(mathrmY) complex random numbers sampled from the standard normal distribution. Next, we only need to note that taking the partial trace of a pure state psirangle is equivalent to calculating AA^dagger where A is a matrix obtained from reshaping psirangle.Sampling random mixed states is implemented using the HilbertSchmidtStates{β, K} type. The meaning of the type parameters is the same as in the Wishart matrices case. We provide additional constructors which set the default values of the parameters HilbertSchmidtStates{β}(d::Int) where β = HilbertSchmidtStates{β, 1}(d), HilbertSchmidtStates(d::Int) = HilbertSchmidtStates{2, 1}(d). The latter one is the most frequent use case. Here is an exampleh = HilbertSchmidtStates(3)\n\nρ = rand(h)\n\ntr(ρ)\n\neigvals(ρ)"
},

{
    "location": "man/random/#Random-quantum-channels-1",
    "page": "Random quantum objects",
    "title": "Random quantum channels",
    "category": "section",
    "text": "Quantum channels are a special subclass of quantum states with constraints imposed on their partial trace as well as trace. Formally, we start with a Ginibre matrix G in mathrmL (mathrmX otimes mathrmY mathrmZ). We obtain a random Choi-Jamiołkowski matrix J_Phi corresponding to a channel Phi as J_Phi = left( mathbbI_mathrmX otimes (mathrmTr_mathrmX GG^dagger)^-12 right) GG^dagger left( mathbbI_mathrmX otimes (mathrmTr_mathrmX GG^dagger)^-12 right)When dim(mathrmZ)=dim(mathrmX) dim(mathrmY) this is known to generate a uniform distribution over the set of quantum channels.The implementation uses the type ChoiJamiolkowskiMatrices{β, K}. The parameters beta and K have the same meaning as in the Wishart matrix case. Additionally here, the constructor ChoiJamiolkowskiMatrices{β, K}(idim::Int, odim::Int)  where {β, K} takes two parameters–the input and output dimension of the channel. As in the previous cases we provide some additional constructors for conveniencefunction ChoiJamiolkowskiMatrices{β}(idim::Int, odim::Int) where β     ChoiJamiolkowskiMatrices{β, 1}(idim, odim) end,function ChoiJamiolkowskiMatrices{β}(d::Int) where β     ChoiJamiolkowskiMatrices{β}(d, d) end,function ChoiJamiolkowskiMatrices(idim::Int, odim::Int)     ChoiJamiolkowskiMatrices{2}(idim, odim) end,function ChoiJamiolkowskiMatrices(d::Int)     ChoiJamiolkowskiMatrices(d, d) end.Here is an example of usagec = ChoiJamiolkowskiMatrices(2, 3)\n\nΦ = rand(c)\n\nptrace(Φ.matrix, [3, 2],[1])Note that the resulting sample is of type DynamicalMatrix."
},

{
    "location": "man/random/#refs_rand-1",
    "page": "Random quantum objects",
    "title": "References",
    "category": "section",
    "text": "[1] B. Collins, I. Nechita, Random matrix techniques in quantum information theory, Journal of Mathematical Physics, 2016;57(1):015215.[2] W. K. Wootters, Random quantum statesy, Foundations of Physics, 1990;20(11):1365–1378.[3] K. Życzkowski, H. J. Sommers, Induced measures in the space of mixed quantum states, Journal of Physics A: Mathematical and General, 2001;34(35):7111.[4] H. J. Sommers, K. Życzkowski, Statistical properties of random density matrices, Journal of Physics A: Mathematical and General, 2004;37(35):8457.[5] Z. Puchała, Ł. Pawela, K. Życzkowski, Distinguishability of generic quantum states, Physical Review A, 2016;93(6):062112.[6] L. Zhang, U. Singh, A. K. Pati, Average subentropy, coherence and entanglement of random mixed   quantum states, Annals of Physics, 2017;377:125–146.[7] L. Zhang, Average coherence and its typicality for random mixed quantum states, Journal of Physics A: Mathematical and Theoretical,   2017;50(15):155303.[8] W. Bruzda, V. Cappellini, H. J. Sommers, K. Życzkowski, Random quantum operations, Physics Letters A,   2009;373(3):320–324.[9] I. Nechita, Z. Puchała, L. Pawela, K. Życzkowski, Almost all quantum channels are equidistant, Journal of Mathematical Physics,   2018;59(5):052201.[10] L. Zhang, J. Wang, Z. Chen, Spectral density of mixtures of random density matrices for qubits, Physics Letters A,   2018;382(23):1516–1523.[11] J. Ginibre, Statistical ensembles of complex, quaternion, and real matrices, Journal of Mathematical Physics, 1965;6(3):440–449.[12] M. L. Mehta, Random matrices. vol. 142., Elsevier; 2004.[13] K. Życzkowski, M. Kuś, Random unitary matrices, Journal of Physics A: Mathematical and General, 1994;27(12):4235.[14] C. Jarlskog, A recursive parametrization of unitary matrices, Journal of Mathematical Physics, 2005;46(10):103508."
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
    "location": "lib/QuantumInformation/#QuantumInformation.SuperOperator-Union{Tuple{T}, Tuple{Function,Int64,Int64}} where T<:(AbstractArray{#s14,2} where #s14<:Number)",
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
    "location": "lib/QuantumInformation/#QuantumInformation.applychannel-Tuple{KrausOperators{#s14} where #s14<:(AbstractArray{#s13,2} where #s13<:Number),AbstractArray{#s12,2} where #s12<:Number}",
    "page": "Library",
    "title": "QuantumInformation.applychannel",
    "category": "method",
    "text": "applychannel(Φ, ρ)\n\n\nΦ: list of vectors.\nρ: input matrix.\n\nReturn application of channel Φonρ`. Kraus representation of quantum channel Phi is a set K_i_iin I of bounded operators on mathcalH such that sum_iin I K_i^dagger K_i = mathcal1. Then Phi(rho)=sum_iin I K_i rho K_i^dagger.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.applychannel-Tuple{Stinespring{#s14} where #s14<:(AbstractArray{#s13,2} where #s13<:Number),AbstractArray{#s12,2} where #s12<:Number}",
    "page": "Library",
    "title": "QuantumInformation.applychannel",
    "category": "method",
    "text": "applychannel(Φ, ρ)\n\n\nΦ: Stinespring representation of quantum channel.\nρ: quantum state.\ndims: dimensions of registers of ρ.\n\nApplication of Stinespring representation of quantum channel into state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.applychannel-Tuple{SuperOperator{#s14} where #s14<:(AbstractArray{#s13,2} where #s13<:Number),AbstractArray{#s12,2} where #s12<:Number}",
    "page": "Library",
    "title": "QuantumInformation.applychannel",
    "category": "method",
    "text": "applychannel(Φ, ρ)\n\n\nΦ: super-operator matrix.\nρ: quantum state.\n\nApplication of super-operator matrix into state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.applychannel-Union{Tuple{T}, Tuple{DynamicalMatrix{#s16} where #s16<:(AbstractArray{#s14,2} where #s14<:Number),AbstractArray{T,2}}} where T<:Number",
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
    "location": "lib/QuantumInformation/#QuantumInformation.bures_angle-Tuple{AbstractArray{#s22,2} where #s22<:Number,AbstractArray{#s21,2} where #s21<:Number}",
    "page": "Library",
    "title": "QuantumInformation.bures_angle",
    "category": "method",
    "text": "bures_angle(ρ, σ)\n\n\nρ: quantum state.\nσ: quantum state.\n\nReturn Bures angle between quantum states ρ and σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.bures_distance-Tuple{AbstractArray{#s22,2} where #s22<:Number,AbstractArray{#s21,2} where #s21<:Number}",
    "page": "Library",
    "title": "QuantumInformation.bures_distance",
    "category": "method",
    "text": "bures_distance(ρ, σ)\n\n\nρ: quantum state.\nσ: quantum state.\n\nReturn Bures distance between quantum states ρ and σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.concurrence-Tuple{AbstractArray{#s23,2} where #s23<:Number}",
    "page": "Library",
    "title": "QuantumInformation.concurrence",
    "category": "method",
    "text": "concurrence(ρ)\n\n\nρ: quantum state.\n\nCalculates the concurrence of a two-qubit system ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.diamond_distance-Union{Tuple{T}, Tuple{DynamicalMatrix{T},DynamicalMatrix{T},Vararg{Any,N} where N}} where T<:(AbstractArray{#s23,2} where #s23<:Number)",
    "page": "Library",
    "title": "QuantumInformation.diamond_distance",
    "category": "method",
    "text": "diamond_distance(Φ1, Φ2, args)\n\n\nΦ1: DynamicalMatrix\nΦ2: DynamicalMatrix\n\nReturn diamond distance between dynamical matrices Φ1 and Φ2.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.fidelity-Tuple{AbstractArray{#s22,2} where #s22<:Number,AbstractArray{#s21,2} where #s21<:Number}",
    "page": "Library",
    "title": "QuantumInformation.fidelity",
    "category": "method",
    "text": "fidelity(ρ, σ)\n\n\nρ: matrix.\nσ: matrix.\n\nReturn fidelity between matrices ρ and σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.fidelity_sqrt-Tuple{AbstractArray{#s22,2} where #s22<:Number,AbstractArray{#s21,2} where #s21<:Number}",
    "page": "Library",
    "title": "QuantumInformation.fidelity_sqrt",
    "category": "method",
    "text": "fidelity_sqrt(ρ, σ)\n\n\nρ: matrix.\nσ: matrix.\n\nReturn square root of fidelity between matrices ρ and σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.gate_fidelity-Tuple{AbstractArray{#s22,2} where #s22<:Number,AbstractArray{#s21,2} where #s21<:Number}",
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
    "location": "lib/QuantumInformation/#QuantumInformation.hs_distance-Tuple{AbstractArray{#s22,2} where #s22<:Number,AbstractArray{#s21,2} where #s21<:Number}",
    "page": "Library",
    "title": "QuantumInformation.hs_distance",
    "category": "method",
    "text": "hs_distance(A, B)\n\n\nA: matrix.\nB: matrix.\n\nReturn Hilbert–Schmidt distance between matrices A and B.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.iscptp-Tuple{KrausOperators{#s22} where #s22<:(AbstractArray{#s23,2} where #s23<:Number)}",
    "page": "Library",
    "title": "QuantumInformation.iscptp",
    "category": "method",
    "text": "iscptp(Φ; atol)\n\n\nΦ: list of Kraus operators.\natol: tolerance of approximation.\n\nChecks if set of Kraus operators fulfill completness relation.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.js_divergence-Tuple{AbstractArray{#s22,2} where #s22<:Number,AbstractArray{#s21,2} where #s21<:Number}",
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
    "location": "lib/QuantumInformation/#QuantumInformation.kl_divergence-Tuple{AbstractArray{#s22,2} where #s22<:Number,AbstractArray{#s21,2} where #s21<:Number}",
    "page": "Library",
    "title": "QuantumInformation.kl_divergence",
    "category": "method",
    "text": "kl_divergence(ρ, σ)\n\n\nρ: quantum state.\nσ: quantum state.\n\nReturn Kullback–Leibler divergence of quantum state ρ with respect to σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.log_negativity-Tuple{AbstractArray{#s23,2} where #s23<:Number,Array{Int64,1},Int64}",
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
    "location": "lib/QuantumInformation/#QuantumInformation.negativity-Tuple{AbstractArray{#s23,2} where #s23<:Number,Array{Int64,1},Int64}",
    "page": "Library",
    "title": "QuantumInformation.negativity",
    "category": "method",
    "text": "negativity(ρ, dims, sys)\n\n\nρ: quantum state.\ndims: dimensions of subsystems.\nsys: transposed subsystem.\n\nReturn negativity of quantum state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.norm_diamond-Union{Tuple{DynamicalMatrix{T}}, Tuple{T}, Tuple{DynamicalMatrix{T},Any}, Tuple{DynamicalMatrix{T},Any,Any}} where T<:(AbstractArray{#s21,2} where #s21<:Number)",
    "page": "Library",
    "title": "QuantumInformation.norm_diamond",
    "category": "method",
    "text": "norm_diamond(Φ)\nnorm_diamond(Φ, method)\nnorm_diamond(Φ, method, eps)\n\n\nΦ: DynamicalMatrix\n\nReturn diamond norm of dynamical matrix Φ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.norm_hs-Tuple{AbstractArray{#s23,2} where #s23<:Number}",
    "page": "Library",
    "title": "QuantumInformation.norm_hs",
    "category": "method",
    "text": "norm_hs(A)\n\n\nA: matrix.\n\nReturn Hilbert–Schmidt norm of matrix A.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.norm_trace-Tuple{AbstractArray{#s23,2} where #s23<:Number}",
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
    "location": "lib/QuantumInformation/#QuantumInformation.ppt-Tuple{AbstractArray{#s23,2} where #s23<:Number,Array{Int64,1},Int64}",
    "page": "Library",
    "title": "QuantumInformation.ppt",
    "category": "method",
    "text": "ppt(ρ, dims, sys)\n\n\nρ: quantum state.\ndims: dimensions of subsystems.\nsys: transposed subsystem.\n\nReturn minimum eigenvalue of positive partial transposition of quantum state ρ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.proj-Tuple{AbstractArray{#s13,1} where #s13<:Number}",
    "page": "Library",
    "title": "QuantumInformation.proj",
    "category": "method",
    "text": "proj(ψ)\n\n\nket: input column vector.\n\nReturn outer product ketranglelangle ket of ket.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ptrace-Tuple{AbstractArray{#s22,1} where #s22<:Number,Array{Int64,1},Int64}",
    "page": "Library",
    "title": "QuantumInformation.ptrace",
    "category": "method",
    "text": "ptrace(ψ, idims, sys)\n\n\nψ: quantum state pure state (ket).\nidims: dimensins of subsystems - only bipartite states accepted.\nsys: traced subsystem.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ptrace-Tuple{AbstractArray{#s23,2} where #s23<:Number,Array{Int64,1},Array{Int64,1}}",
    "page": "Library",
    "title": "QuantumInformation.ptrace",
    "category": "method",
    "text": "ptrace(ρ, idims, isystems)\n\n\nρ: quantum state.\nidims: dimensins of subsystems.\nisystems: traced subsystems.\n\nReturn partial trace of matrix ρ over the subsystems determined by isystems.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ptrace-Tuple{AbstractArray{#s23,2} where #s23<:Number,Array{Int64,1},Int64}",
    "page": "Library",
    "title": "QuantumInformation.ptrace",
    "category": "method",
    "text": "ptrace(ρ, idims, sys)\n\n\nρ: quantum state.\nidims: dimensins of subsystems.\nsys: traced subsystem.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ptranspose-Tuple{AbstractArray{#s22,2} where #s22<:Number,Array{Int64,1},Array{Int64,1}}",
    "page": "Library",
    "title": "QuantumInformation.ptranspose",
    "category": "method",
    "text": "ptranspose(ρ, idims, isystems)\n\n\nρ: quantum state.\nidims: dimensins of subsystems.\nisystems: transposed subsystems.\n\nReturn partial transposition of matrix ρ over the subsystems determined by isystems.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.ptranspose-Tuple{AbstractArray{#s23,2} where #s23<:Number,Array{Int64,1},Int64}",
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
    "location": "lib/QuantumInformation/#QuantumInformation.relative_entropy-Tuple{AbstractArray{#s22,2} where #s22<:Number,AbstractArray{#s21,2} where #s21<:Number}",
    "page": "Library",
    "title": "QuantumInformation.relative_entropy",
    "category": "method",
    "text": "relative_entropy(ρ, σ)\n\n\nρ: quantum state.\nσ: quantum state.\n\nReturn quantum relative entropy of quantum state ρ with respect to σ.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.res-Tuple{AbstractArray{#s13,2} where #s13<:Number}",
    "page": "Library",
    "title": "QuantumInformation.res",
    "category": "method",
    "text": "res(ρ)\n\n\nρ: input matrix.\n\nReturns vec(ρ.T). Reshaping maps     matrix ρ into a vector row by row.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.reshuffle-Tuple{AbstractArray{#s22,2} where #s22<:Number}",
    "page": "Library",
    "title": "QuantumInformation.reshuffle",
    "category": "method",
    "text": "reshuffle(ρ)\n\n\nρ: reshuffled matrix.\n\nPerforms reshuffling of indices of a matrix.   Given multiindexed matrix M_(mμ)(nν) it returns   matrix M_(mn)(μν).\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.shannon_entropy-Tuple{AbstractArray{#s23,1} where #s23<:Real}",
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
    "location": "lib/QuantumInformation/#QuantumInformation.superfidelity-Tuple{AbstractArray{#s22,2} where #s22<:Number,AbstractArray{#s21,2} where #s21<:Number}",
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
    "location": "lib/QuantumInformation/#QuantumInformation.unres-Tuple{AbstractArray{#s13,1} where #s13<:Number}",
    "page": "Library",
    "title": "QuantumInformation.unres",
    "category": "method",
    "text": "unres(ρ)\n\n\nϕ: input matrix.\n\nReturn de-reshaping of the vector into a matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.vonneumann_entropy-Tuple{LinearAlgebra.Hermitian{#s23,S} where S<:(AbstractArray{#s576,2} where #s576<:#s23) where #s23<:Number}",
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
    "text": "werner_state(d, α)\n\n\nd: length of the vector.\nα: real number from [0, 1].\n\nReturns Werner state given by fracalphadleft(sum_i=0^sqrtd-1iirangleright) left(sum_i=0^sqrtd-1langle iiright)+ frac1-alphadsum_i=0^d-1iranglelangle i.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{DynamicalMatrix{T1}},KrausOperators{T2}}} where T2<:(AbstractArray{#s14,2} where #s14<:Number) where T1<:(AbstractArray{#s16,2} where #s16<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\n?: type.\nΦ: list of Kraus operators.\n\nTransforms list of Kraus operators into dynamical matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{DynamicalMatrix{T1}},SuperOperator{T2}}} where T2<:(AbstractArray{#s14,2} where #s14<:Number) where T1<:(AbstractArray{#s16,2} where #s16<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\n?: type.\nΦ: super-operator matrix.\n\nTransforms super-operator matrix into dynamical matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{KrausOperators{T1}},DynamicalMatrix{T2}}} where T2<:(AbstractArray{#s12,2} where #s12<:Number) where T1<:(AbstractArray{#s13,2} where #s13<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\n?: type.\nΦ: dynamical matrix.\n\nTransforms dynamical matrix into list of Kraus operators.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{KrausOperators{T1}},SuperOperator{T2}}} where T2<:(AbstractArray{#s14,2} where #s14<:Number) where T1<:(AbstractArray{#s16,2} where #s16<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\n?: type.\nΦ: super-operator matrix.\n\nTransforms super-operator matrix into list of Kraus operators.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{Stinespring{T1}},DynamicalMatrix{T2}}} where T2<:(AbstractArray{#s14,2} where #s14<:Number) where T1<:(AbstractArray{#s16,2} where #s16<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\n?: type.\nΦ: dynamical matrix.\n\nTransforms dynamical matrix into Stinespring representation of quantum channel.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{Stinespring{T1}},KrausOperators{T2}}} where T2<:(AbstractArray{#s12,2} where #s12<:Number) where T1<:(AbstractArray{#s13,2} where #s13<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\n?: type.\nΦ: list of Kraus operators.\n\nTransforms list of Kraus operators into Stinespring representation of quantum channel.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{Stinespring{T1}},SuperOperator{T2}}} where T2<:(AbstractArray{#s14,2} where #s14<:Number) where T1<:(AbstractArray{#s16,2} where #s16<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\n?: type.\nΦ: super-operator matrix.\n\nTransforms super-operator matrix into Stinespring representation of quantum channel.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{SuperOperator{T1}},DynamicalMatrix{T2}}} where T2<:(AbstractArray{#s14,2} where #s14<:Number) where T1<:(AbstractArray{#s16,2} where #s16<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\n?: type.\nΦ: dynamical matrix.\n\nTransforms dynamical matrix into super-operator matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#Base.convert-Union{Tuple{T2}, Tuple{T1}, Tuple{Type{SuperOperator{T1}},KrausOperators{T2}}} where T2<:(AbstractArray{#s14,2} where #s14<:Number) where T1<:(AbstractArray{#s16,2} where #s16<:Number)",
    "page": "Library",
    "title": "Base.convert",
    "category": "method",
    "text": "convert(?, Φ)\n\n\n?: type.\nΦ: list of Kraus operators.\n\nTransforms list of Kraus operators into super-operator matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/QuantumInformation/#QuantumInformation.purity-Tuple{AbstractArray{#s23,2} where #s23<:Number}",
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

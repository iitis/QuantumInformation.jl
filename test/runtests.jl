using QI
if VERSION<=v"0.7"
    const ComplexF64 = Complex128
    using Base.Test
else
    using LinearAlgebra
    using Test
end

my_tests = ["utils.jl", "base.jl", "ptrace.jl", "ptranspose.jl", "reshuffle.jl",
"randommatrix.jl", "randomstate.jl", "functionals.jl",
"gates.jl"]

for my_test in my_tests
    include(my_test)
end

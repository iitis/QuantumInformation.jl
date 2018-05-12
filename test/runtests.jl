using QI
if VERSION<=v"0.7"
    const ComplexF64 = Complex128
    using Base.Test
else
    using LinearAlgebra
    using Test
end

my_tests = ["utils.jl", "base.jl", "randommatrix.jl", "randomstate.jl"]

for my_test in my_tests
    include(my_test)
end

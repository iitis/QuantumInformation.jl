using QI
if VERSION<=v"0.7"
    const ComplexF64 = Complex128
    using Base.Test
else
    using LinearAlgebra
    using Test
end

my_tests = ["utils.jl", "base.jl", "randommatrix.jl", "randomstate.jl"]

println("Running tests:")

for my_test in my_tests
    println(" * $(my_test)")
    include(my_test)
end

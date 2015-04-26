using QI
⊗ = kron
using Base.Test

my_tests = ["utils.jl", "base.jl", "randommatrix.jl", "randomstate.jl"]

println("Running tests:")

for my_test in my_tests
    println(" * $(my_test)")
    include(my_test)
end
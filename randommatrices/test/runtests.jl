using Test
using Random

include("../src/RandomMatrices.jl")
using .RandomMatrices

my_tests = ["circular.jl", "ginibre.jl", "wigner.jl", "wishart.jl"]

for my_test in my_tests
    include(my_test)
end

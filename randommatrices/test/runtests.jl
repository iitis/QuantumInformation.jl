using Test
using Random
using StatsBase

include("../src/RandomMatrices.jl")
using .RandomMatrices

my_tests = ["circular.jl", "ginibre.jl", "wigner.jl", "wishart.jl"]

for my_test in my_tests
    include(my_test)
end

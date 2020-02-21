using Test
using Random
using StatsBase
using LinearAlgebra
using CuArrays

include("../../randommatrices/src/RandomMatrices.jl")
using ..RandomMatrices

include("../src/CuRandomMatrices.jl")
using .CuRandomMatrices

my_tests = ["circular.jl", "ginibre.jl", "wigner.jl", "wishart.jl"]

for my_test in my_tests
    include(my_test)
end

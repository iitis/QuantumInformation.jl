using Test
using LinearAlgebra

include("../src/ElementaryArrays.jl")
using .ElementaryArrays

my_tests = ["elementaryarrays.jl"]

for my_test in my_tests
    include(my_test)
end

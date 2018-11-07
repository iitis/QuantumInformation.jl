using QI
using Random

using LinearAlgebra
using SparseArrays
using Test

my_tests = ["utils.jl", "base.jl", "ptrace.jl", "ptranspose.jl", "reshuffle.jl",
            "channels.jl", "functionals.jl", "gates.jl", "permute_systems.jl",
            "randomqobjects.jl"]

push!(my_tests, "convex.jl") # Convex.jl does not support julia 0.7 yet

println(my_tests)

for my_test in my_tests
    include(my_test)
end

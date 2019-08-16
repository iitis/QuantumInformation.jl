using QuantumInformation
using Random

using LinearAlgebra
using SparseArrays
using Test

my_tests = ["utils.jl", "base.jl", "ptrace.jl", "ptranspose.jl", "reshuffle.jl",
            "channels.jl", "functionals.jl", "gates.jl", "matrixbases.jl",
            "permute_systems.jl", "randomqobjects.jl", "convex.jl"]

for my_test in my_tests
    include(my_test)
end

using QuantumInformation
using Random

using LinearAlgebra
# using SparseArrays
using Test
using Aqua
using JET

my_tests = [
    "utils.jl",
    "base.jl",
    "ptrace.jl",
    "ptranspose.jl",
    "reshuffle.jl",
    "channels.jl",
    "functionals.jl",
    "gates.jl",
    "matrixbases.jl",
    "permute_systems.jl",
    "randomqobjects.jl",
    "convex.jl"
    ]
for my_test in my_tests
    include(my_test)
end

@testset "Aqua.jl" begin
    Aqua.test_all(QuantumInformation)
end

@testset "JET.jl" begin
    JET.test_package(QuantumInformation; target_defined_modules=true)
end

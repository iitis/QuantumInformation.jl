using QuantumInformation
using Random

using LinearAlgebra
# using SparseArrays
using Test
using Aqua
using JET

@testset verbose=true "QuantumInformation.jl" begin
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
        @testset "$my_test" begin
            include(my_test)
        end
    end

    @testset "Aqua.jl" begin
        Aqua.test_all(QuantumInformation)
    end

    @testset "JET.jl" begin
        JET.test_package(QuantumInformation; target_defined_modules=true)
    end
end

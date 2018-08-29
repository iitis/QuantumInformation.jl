using QI
using Random

if VERSION<v"0.7.0"
    const ComplexF64 = Complex128
    using Base.Test
else
    using LinearAlgebra
    using SparseArrays
    using Test
end

my_tests = ["utils.jl", "base.jl", "ptrace.jl", "ptranspose.jl", "reshuffle.jl",
            "channels.jl", "functionals.jl", "gates.jl", "permute_systems.jl",
            "randomqobjects.jl"]

if VERSION<v"0.7.0"
    push!(my_tests, "convex.jl")
else
    # Convex.jl does not support julia 0.7 yet
end

for my_test in my_tests
    include(my_test)
end

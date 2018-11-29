using QuantumInformation
using LinearAlgebra
using Statistics
using JLD2
using FileIO
using ArgParse

function comptime(ccalc::Function, steps::Int, dims::Vector)
  ccalc(steps, 2)
  comptimes = zeros(length(dims))
  for (i, d) in enumerate(dims)
    tic = time_ns()
    ccalc(steps, d)
    toc = time_ns()
    comptimes[i] = Float64(toc-tic)/1.0e9
  end
  comptimes
end

function random_unitary(steps::Int, d::Int)
  dist = CUE(d)
  for i=1:steps
    U = rand(dist)
  end
end

function random_pure_state(steps::Int, d::Int)
  dist = HaarKet(d)
  for i=1:steps
    ψ = rand(dist)
  end
end

function random_mixed_state(steps::Int, d::Int)
  dist = HilbertSchmidtStates(d)
  for i=1:steps
    ρ =rand(dist)
  end
end

function random_channel(steps::Int, d::Int)
  dist = ChoiJamiolkowskiMatrices(round(Int, sqrt(d)))
  for i=1:steps
    Φ = convert(SuperOperator{Matrix{ComplexF64}}, rand(dist))
  end
end

function trace_distance_random(steps::Int, d::Int)
  dist = HilbertSchmidtStates(d)
  for i=1:steps
    trace_distance(rand(dist), rand(dist))
  end
end

function trace_distance_max_mixed(steps::Int, d::Int)
  dist = HilbertSchmidtStates(d)
  ρ = 𝕀(d)/d
  for i=1:steps
    trace_distance(rand(dist), ρ)
  end
end

function entropy_stationary(steps::Int, d::Int)
    dist = ChoiJamiolkowskiMatrices(round(Int, sqrt(d)))
    for i=1:steps
      Φ = rand(dist)
      F = eigen(convert(SuperOperator{Matrix{ComplexF64}}, Φ).matrix)
      idx = findfirst(x->isapprox(x, 1), F.values)
      ρ = unres(F.vectors[:, idx])
      ρ = Hermitian(ρ ./ tr(ρ))
      vonneumann_entropy(ρ)
    end
end

function savect(steps::Int, dims::Vector)
  filename = replace("res/$(steps)_$(dims)_random.jld2", "["=>"")
  filename = replace(filename, "]"=>"")
  filename = replace(filename, " "=>"")
  cases = [(random_unitary, "random_unitary"), (random_pure_state, "random_pure_state"),
  (random_mixed_state, "random_mixed_state"), (random_channel, "random_channel"),
  (entropy_stationary, "entropy_stationary"), (trace_distance_max_mixed, "trace_distance_max_mixed"),
  (trace_distance_random, "trace_distance_random")]
  compt = Dict{String, Any}("steps" => steps, "dims" => dims)
  for (ccalc, label) in cases
    println(label)
    push!(compt, label => comptime(ccalc, steps, dims) ./ steps)
  end
  println(compt)
  save(filename, compt)
end

function main(args)
  s = ArgParseSettings("description")
  @add_arg_table s begin
      "--steps", "-s"
        help = "number of samples per dimension"
        default = 10
        arg_type = Int
      "--dims", "-d"
        help = "dimensions"
        nargs = '*'
        default = [4, 16, 64]
        arg_type = Int
    end
  parsed_args = parse_args(s)
  steps = parsed_args["steps"]
  dims = parsed_args["dims"]
  savect(steps, dims)
end

main(ARGS)

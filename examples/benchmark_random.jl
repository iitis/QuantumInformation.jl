using QI
using LinearAlgebra
using Statistics
using JLD2
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
    ρ = rand(dist)
  end
end

function random_channel(steps::Int, d::Int)
  dist = ChoiJamiolkowskiMatrices(round(Int, sqrt(d)))
  for i=1:steps
    Φ = rand(dist)
  end
end

function savect(steps::Int, dims::Vector)
  filename = replace("res/$(steps)_$(dims)_random.jld2", "["=>"")
  filename = replace(filename, "]"=>"")
  cases = [(random_unitary, "random_unitary"), (random_pure_state, "random_pure_state"), (random_mixed_state, "random_mixed_state"), (random_channel, "random_channel")]
  compt = Dict{String, Any}("steps" => steps, "dims" => dims)
  for (ccalc, label) in cases
    println(label)
    push!(compt, label => comptime(ccalc, steps, dims) ./ steps)
  end
  println(compt)
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

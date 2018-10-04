using JLD2
using FileIO
using NPZ
using Plots
##
filename_base = "examples/res/100_4,16,64,256_random"
##
data_qutip = NPZ.npzread(filename_base*".npz")
data_julia = load(filename_base*".jld2")
##
labels = keys(data_julia)
results = Dict()
p = undef
for l in labels
    if l == "dims" || l == "steps" continue end
    println(l)
    results[l] = data_qutip[l] ./ data_julia[l]
    # p = plot(data_qutip["dims"], results[l], title=l)
    p=plot(data_qutip["dims"], [data_qutip[l], data_julia[l]], title=l, label=["qutip", "julia"])
    display(p)
end
# l = collect(labels)[3]
# plot(data_qutip["dims"], results[l])

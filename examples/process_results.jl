using JLD2
using FileIO
using NPZ
using PyCall
@pyimport matplotlib as mpl
mpl.rc("text", usetex=true)
@pyimport pylab as pl
using LaTeXStrings
##
filename_base = "res/100_4,16,64,256,1024_random"
##
data_qutip = NPZ.npzread(filename_base*".npz")
data_julia = load(filename_base*".jld2")
##
labels = keys(data_julia)
results = Dict()
for l in labels
    if l == "dims" || l == "steps" continue end
    println(l)
    results[l] = data_qutip[l] ./ data_julia[l]

    fig, ax = pl.subplots(figsize=(4,3))
    ax[:semilogy](data_qutip["dims"], data_qutip[l], "-o",
    label="\$\\mathrm{QuTiP}\$")
    ax[:semilogy](data_qutip["dims"], data_julia[l], "-^",
    label="\$\\mathrm{QuantumInformation.jl}\$")
    fig[:legend](bbox_to_anchor = (0.9,0.6))
    ax[:set_xticks]([4, 64, 256, 1024])
    ax[:set_xlabel]("\$d\$")
    ax[:set_ylabel]("\$t\\ (s)\$")
    pl.savefig("$(filename_base)_$(l).pdf", bbox_inches="tight")
    pl.close(fig)
end

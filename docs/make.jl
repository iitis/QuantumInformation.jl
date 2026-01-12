using Documenter, QuantumInformation

format = Documenter.HTML(
    edit_link="master",
    prettyurls=get(ENV, "CI", nothing) == "true",
    assets=["assets/favicon.ico"],
)

makedocs(
    clean=true,
    format=format,
    sitename="QuantumInformation.jl",
    authors="Piotr Gawron, Dariusz Kurzyk, Łukasz Pawela",
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
            "man/quickstart.md",
            "man/vectors.md",
            "man/states.md",
            "man/functionals.md",
            "man/measurement.md",
            "man/channels.md",
            "man/random.md",
        ],
        "Library" => "lib/QuantumInformation.md",
        # Any[
        #     "lib/QuantumInformation.md",
        #     "lib/content/base.md",
        #     "lib/content/gates.md",
        #     "lib/content/randommatrix.md",
        #     "lib/content/randomstate.md",
        #     "lib/content/utils.md"
        # ]
    ],
)

deploydocs(
    target="build",
    repo="github.com/iitis/QuantumInformation.jl.git",
)

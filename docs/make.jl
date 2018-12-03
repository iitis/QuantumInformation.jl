using Documenter, QuantumInformation

makedocs(
    clean = true,
    format = :html,
    sitename = "QuantumInformation.jl",
    authors = "Piotr Gawron, Dariusz Kurzyk, Łukasz Pawela",
    assets = ["assets/favicon.ico"],
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "man/quickstart.md",
            "man/vectors.md",
            "man/states.md",
            "man/functionals.md",
            "man/measurement.md",
            "man/random.md"
        ],
        "Library" => "lib/QuantumInformation.md"
        # Any[
        #     "lib/QuantumInformation.md",
        #     "lib/content/base.md",
        #     "lib/content/gates.md",
        #     "lib/content/randommatrix.md",
        #     "lib/content/randomstate.md",
        #     "lib/content/utils.md"
        # ]
    ]
)

deploydocs(
    deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    target = "build",
    repo = "github.com/ZKSI/QuantumInformation.jl.git"
)

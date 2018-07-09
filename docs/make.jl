using Documenter, QI

# Pusty makedocs generuje dokumentacje w domyslnym
# formacie markdown.
# makedocs(
#
# )
makedocs(
    modules = [QI],
    clean = false,
    format = :html,
    sitename = "QI.jl",
    authors = "Piotr Gawron, Dariusz Kurzyk, Łukasz Pawela",
    assets = ["assets/favicon.ico"],
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "man/quickstart.md",
            "man/vectors.md",
            "man/states.md",
            "man/functionals.md",
            "man/random.md"
        ],
        "Library" => "lib/QI.md"
        # Any[
        #     "lib/QI.md",
        #     "lib/content/base.md",
        #     "lib/content/gates.md",
        #     "lib/content/randommatrix.md",
        #     "lib/content/randomstate.md",
        #     "lib/content/utils.md"
        # ]
    ],
    html_canonical = "https://zksi.github.io/QI.jl/latest/",
)

deploydocs(
    deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    target = "build",
    repo = "github.com/ZKSI/QI.jl.git",
    julia = "0.6.3",
    deps = nothing,
    make = nothing
)

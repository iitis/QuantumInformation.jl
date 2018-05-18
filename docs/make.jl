using Documenter, QI

makedocs()

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/ZKSI/QI.jl.git"
)

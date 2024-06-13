using Documenter
using Cosmology
using Documenter.Remotes: GitHub


DocMeta.setdocmeta!(Cosmology, :DocTestSetup, :(using Cosmology); recursive = true)

include("pages.jl")

makedocs(;
    modules = [Cosmology],
    authors = "Julia Astro",
    repo = GitHub("JuliaAstro/Cosmology.jl"),
    sitename = "Cosmology.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliaastro.github.io/Cosmology.jl",
        assets = String[],
   ),
    pages=pages,
)

deploydocs(;
    repo = "github.com/JuliaAstro/Cosmology.jl",
)

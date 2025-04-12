using Cosmology
using Documenter
using Documenter.Remotes: GitHub
using DocumenterCitations


DocMeta.setdocmeta!(Cosmology, :DocTestSetup, :(using Cosmology); recursive = true)

include("pages.jl")

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

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
    plugins = [bib],
    pages,
)

deploydocs(;
    repo = "github.com/JuliaAstro/Cosmology.jl",
    push_preview = true,
)

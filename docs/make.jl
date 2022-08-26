using Documenter
using Cosmology


DocMeta.setdocmeta!(Cosmology, :DocTestSetup, :(using Cosmology); recursive = true)

include("pages.jl")

makedocs(;
    modules = [Cosmology],
    authors = "Julia Astro",
    repo = "https://github.com/JuliaAstro/Cosmology.jl/blob/{commit}{path}#L{line}",
    sitename = "Cosmology.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliaastro.github.io/Cosmology.jl",
        assets = String[],),
    pages=pages,)

deploydocs(;
    repo = "github.com/JuliaAstro/Cosmology.jl",)

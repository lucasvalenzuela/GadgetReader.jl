using GadgetGalaxies
using Documenter

DocMeta.setdocmeta!(GadgetGalaxies, :DocTestSetup, :(using GadgetGalaxies); recursive=true)

makedocs(;
    modules=[GadgetGalaxies],
    authors="Lucas Valenzuela <lucasvalenzuela@users.noreply.github.com> and contributors",
    repo="https://github.com/lucasvalenzuela/GadgetGalaxies.jl/blob/{commit}{path}#{line}",
    sitename="GadgetGalaxies.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lucasvalenzuela.github.io/GadgetGalaxies.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lucasvalenzuela/GadgetGalaxies.jl",
)

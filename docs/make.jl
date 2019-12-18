using Documenter, Graphene

makedocs(;
    modules=[Graphene],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/chenspc/Graphene.jl/blob/{commit}{path}#L{line}",
    sitename="Graphene.jl",
    authors="Chen Huang",
    assets=String[],
)

deploydocs(;
    repo="github.com/chenspc/Graphene.jl",
)

using Documenter, Graphenet

makedocs(;
    modules=[Graphenet],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/chenspc/Graphenet.jl/blob/{commit}{path}#L{line}",
    sitename="Graphenet.jl",
    authors="Chen Huang",
    assets=String[],
)

deploydocs(;
    repo="github.com/chenspc/Graphenet.jl",
)

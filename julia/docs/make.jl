using Documenter
using KrySBAS

makedocs(
    sitename = "KrySBAS.jl",
    modules  = [KrySBAS],
    authors  = "CC&MA – NIDTec – FP – UNA",
    format   = Documenter.HTML(
        canonical        = "https://nidtec-una.github.io/krysbas-dev/",
        edit_link        = "main",
        assets           = String[],
    ),
    pages = [
        "Home"    => "index.md",
        "Solvers" => "solvers.md",
    ],
    warnonly = [:missing_docs],
)

deploydocs(
    repo      = "github.com/nidtec-una/krysbas-dev.git",
    devbranch = "main",
    push_preview = true,
)

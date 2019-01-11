using Documenter, VarianceComponentModels

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = :html,
    sitename = "VarianceComponentModels.jl",
    modules = [VarianceComponentModels], 
    pages = Any[
        "Home" => "index.md",
        "Manual" => Any[
            "man/mle_reml.md",
            "man/heritability.md"
        ],
        "API" => "man/api.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/VarianceComponentModels.jl.git",
    target = "build",
    osname = "linux",
    julia = "1.0",
    deps = nothing,
    make = nothing)

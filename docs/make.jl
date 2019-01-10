using Documenter, VarianceComponentModels

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = :html,
    sitename = "VarianceComponentModels",
    modules = [VarianceComponentModels]
)

deploydocs(
    repo   = "github.com/OpenMendel/VarianceComponentModels.jl.git",
    target = "build"
)

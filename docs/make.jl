using Documenter, VarianceComponentModels

ENV["DOCUMENTER_DEBUG"] = "true"
makedocs()
deploydocs(
  deps   = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
  repo   = "github.com:OpenMendel/VarianceComponentModels.jl.git",
  julia  = "nightly",
  osname = "osx"
  )

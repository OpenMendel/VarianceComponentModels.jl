module VarianceComponentModels

using MathProgBase
using Ipopt
#using Mosek
using KNITRO

# utilities for multivariate calculus
include("multivariate_calculus.jl")
# source for fitting models with 2 variance components
include("two_variance_component.jl")

end # module

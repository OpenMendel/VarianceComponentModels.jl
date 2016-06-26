# VarianceComponentModels

[![Build Status](https://travis-ci.org/Hua-Zhou/VarianceComponentModels.jl.svg?branch=master)](https://travis-ci.org/Hua-Zhou/VarianceComponentModels.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/huae8ed3g4dxhq73/branch/master?svg=true)](https://ci.appveyor.com/project/Hua-Zhou/variancecomponentmodels-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/Hua-Zhou/VarianceComponentModels.jl/badge.svg?branch=master)](https://coveralls.io/github/Hua-Zhou/VarianceComponentModels.jl?branch=master)
[![codecov](https://codecov.io/gh/Hua-Zhou/VarianceComponentModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Hua-Zhou/VarianceComponentModels.jl)

VarianceComponentModels provides computational routines for fitting and testing variance component models.

The response matrix `Y` is assumed to be multivariate normal with mean `X B` and covariance matrix `Σ1 ⊗ V1 + ⋯ + Σm ⊗ Vm`. Here **data** is  
* `Y`: `n x d` response matrix  
* `X`: `n x p` covariate matrix  
* `V1, ..., Vm`: `m` `n x n` covariance matrices
**Parameters** are  
* `B`: `p x d` parameter matrix  
* `Σ1, ..., Σm`: `m` `d x d` variance components

## Installation

Use the Julia package manager to install VarianceComponentModels.jl:

    Pkg.clone("git@github.com:Hua-Zhou/VarianceComponentModels.jl.git")

## Documentation

View the [tutorial]().

## OpenMendel

VarianceComponentModels is one component of the umbrella [OpenMendel]() project. See [cite]() if you use the code in your research.    

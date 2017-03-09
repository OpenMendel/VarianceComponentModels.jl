# VarianceComponentModels.jl

*Utilities for fitting and testing variance component models*

VarianceComponentModels.jl implements computation routines for fitting and testing variance component model of form

$\text{vec}(Y) \sim \text{Nomral}(X B, \Sigma_1 \otimes V_1 + \cdots + \Sigma_m \otimes V_m),$

where $\otimes$ is the [Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product).

In this model, **data** is represented by  

* `Y`: `n x d` response matrix  
* `X`: `n x p` covariate matrix  
* `V=(V1,...,Vm)`: a tuple `m` `n x n` covariance matrices  

and **parameters** are  

* `B`: `p x d` mean parameter matrix  
* `Σ=(Σ1,...,Σm)`: a tuple of `m` `d x d` variance components  

## Package Features

- Maximum likelihood estimation (MLE) and restricted maximum likelihood estimation (REML) of mean parameters `B` and variance component parameters `Σ`   
- Allow constrains in the mean parameters `B`  
- Choice of optimization algorithms: [Fisher scoring](https://books.google.com/books?id=QYqeYTftPNwC&lpg=PP1&pg=PA142#v=onepage&q&f=false) and [minorization-maximization algorithm](http://arxiv.org/abs/1509.07426)  
- [Heritability Analysis](@ref) in genetics  

## Installation

Use the Julia package manager to install VarianceComponentModels.jl.
```julia
Pkg.clone("https://github.com/OpenMendel/VarianceComponentModels.jl.git")
```
This package supports Julia `0.4`.

## Manual Outline

```@contents
Pages = [
    "man/mle_reml.md",
    "man/heritability.md",
]
Depth = 2
```

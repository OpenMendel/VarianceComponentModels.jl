
# API

Documentation for `VarianceComponentModels.jl`'s types and methods.

## Index

```@index
Pages = ["api.md"]
```

## Types

```@docs
VarianceComponentModel
VarianceComponentVariate
TwoVarCompModelRotate
TwoVarCompVariateRotate
```

## Functions

```@docs
eltype
length
size
nmeanparams
nvarparams
nparams
cov!
VarianceComponentModels.cov
mean!
mean
residual
```

```@docs
mle_fs!
mle_mm!
fit_mle!
fit_reml!
heritability
```

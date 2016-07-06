module VarianceComponentModels

using MathProgBase, Ipopt, KNITRO#, Mosek, Gurobi
import Base: eltype, length, size, mean, mean!, cov
export VarianceComponentModel, VarianceComponentVariate,
  TwoVarCompModelRotate, TwoVarCompVariateRotate, residual,
  nvarcomps, nmeanparams, nvarparams, nparams,
  mean!, mean, cov!, cov

"""
`VarianceComponentModel` stores the model parameters of a variance component
model. `B` is `p x d` mean parameters and `Σ` is a tuple of `d x d` variance
component parameters.
"""
type VarianceComponentModel{T <: AbstractFloat, M,
  BT <: AbstractVecOrMat, ΣT <: AbstractMatrix}
  # model parameters
  B::BT
  Σ::NTuple{M, ΣT}
  # constraints: A * vec(B) sense b, lb ≤ vec(B) ≤ ub
  A::Matrix{T}
  sense::Union{Char, Vector{Char}}
  b::Union{T, Vector{T}}
  lb::Union{T, Vector{T}}
  ub::Union{T, Vector{T}}
  # inner constructor
  function VarianceComponentModel(B::AbstractVecOrMat{T},
    Σ::NTuple{M, AbstractMatrix{T}}, A::Matrix{T}, sense::Union{Char, Vector{Char}},
    b::Union{T, Vector{T}}, lb::Union{T, Vector{T}}, ub::Union{T, Vector{T}})
    new(B, Σ, A, sense, b, lb, ub)
  end
end

"""
Default constructor of `VarianceComponentModel` type.
"""
function VarianceComponentModel{M}(
  B::AbstractVecOrMat,
  Σ::NTuple{M, AbstractMatrix},
  A::Matrix = zeros(eltype(B), 0, length(B)),
  sense::Union{Char, Vector{Char}} = Array{Char}(0),
  b::Union{Number, Vector} = zeros(eltype(B), 0),
  lb::Union{Number, Vector} = convert(eltype(B), -Inf),
  ub::Union{Number, Vector} = convert(eltype(B), Inf),
  )

  VarianceComponentModel{eltype(B), M, typeof(B), eltype(Σ)}(B, Σ, A, sense, b,
    lb, ub)
end

"""
Construct a `VarianceComponentModel` instance from `Σ` alone. `B` is treated empty.
"""
function VarianceComponentModel{M}(Σ::NTuple{M, AbstractMatrix})
  B = zeros(eltype(Σ[1]), 0, size(Σ[1], 1))
  VarianceComponentModel(B, Σ)
end

"""
`TwoVarCompModelRotate` stores the rotated two variance component model. If
`B, Σ=(Σ1, Σ2)` is a two variane component model, then
`(eigval, eigvec) = eig(Σ1, Σ2)`, `Brot = B * eigvec`, and `logdetΣ2 = logdet(Σ2)`.
"""
immutable TwoVarCompModelRotate{T <: AbstractFloat, BT <: AbstractVecOrMat}
  Brot::BT
  eigval::Vector{T}
  eigvec::Matrix{T}
  logdetΣ2::T
  function TwoVarCompModelRotate(Brot::AbstractVecOrMat{T},
    eigval::Vector{T}, eigvec::Matrix{T}, logdetΣ2::T)
    new(Brot, eigval, eigvec, logdetΣ2)
  end
end

"""
Default constructor of `TwoVarCompModelRotate` type.
"""
function TwoVarCompModelRotate(
  Brot::AbstractVecOrMat,
  eigval::Vector,
  eigvec::Matrix,
  logdetΣ2
  )

  TwoVarCompModelRotate{eltype(Brot), typeof(Brot)}(Brot, eigval, eigvec, logdetΣ2)
end

"""
Constructor of a `TwoVarCompModelRotate` instance from a
`VarianceComponentModel` instance.
"""
function TwoVarCompModelRotate{T, BT, ΣT}(
  vcm::VarianceComponentModel{T, 2, BT, ΣT}
  )

  # generalized eigenvalue decomposition of (Σ1, Σ2)
  F = eigfact(Symmetric(vcm.Σ[1]), Symmetric(vcm.Σ[2])) # output is of BlasFloat type
  λ = convert(Vector{T}, F.values)
  Φ = convert(Matrix{T}, F.vectors)
  # correct negative eigenvalues due to roundoff
  map!(x -> max(x, zero(T)), λ)
  Brot = isempty(vcm.B) ? Array{T}(size(vcm.B)) : vcm.B * Φ
  logdetΣ2 = convert(T, logdet(vcm.Σ[2]))
  TwoVarCompModelRotate{T, BT}(Brot, λ, Φ, logdetΣ2)
end

"""
`VarianceComponentVariate` stores the data of a variance component
model. `Y` is `n x d` responses, `X` is `n x p` predictors, and `V` is a tuple
of `n x n` covariance matrices.
"""
immutable VarianceComponentVariate{T <: AbstractFloat, M,
  YT <: AbstractVecOrMat, XT <: AbstractVecOrMat, VT <: AbstractMatrix}
  # data
  Y::YT
  X::XT
  V::NTuple{M, VT}
  # inner constructor
  function VarianceComponentVariate(Y::AbstractVecOrMat{T},
    X::AbstractVecOrMat{T}, V::NTuple{M, AbstractMatrix{T}})
    new(Y, X, V)
  end
end

"""
Default constructor of `VarianceComponentVariate` type.
"""
function VarianceComponentVariate{M}(
  Y::AbstractVecOrMat,
  X::AbstractVecOrMat,
  V::NTuple{M, AbstractMatrix}
  )

  VarianceComponentVariate{eltype(Y), M, typeof(Y), typeof(X), eltype(V)}(Y, X, V)
end

"""
Default constructor of a `VarianceComponentVariate` instance from `Y` and `V`
alone. `X` is created empty.
"""
function VarianceComponentVariate{M}(
  Y::AbstractVecOrMat,
  V::NTuple{M, AbstractMatrix}
  )

  X = zeros(eltype(Y), size(Y, 1), 0)
  VarianceComponentVariate{eltype(Y), M, typeof(Y), typeof(X), eltype(V)}(Y, X, V)
end

"""
`TwoVarCompVariateRotate` stores the rotated two variance component data. If
`Y, X, V=(V1, V2)` is the two variane component data, then
`(eigval, eigvec) = eig(V1, V2)`, `Yrot = eigvec' * Y`, `Xrot = eigvec' * X`,
and `logdetV2 = logdet(V2)`.

"""
immutable TwoVarCompVariateRotate{T <: AbstractFloat, YT <: AbstractVecOrMat,
  XT <: AbstractVecOrMat}
  # data
  Yrot::YT
  Xrot::XT
  eigval::Vector{T}
  logdetV2::T
  # inner constructor
  function TwoVarCompVariateRotate(Yrot::AbstractVecOrMat{T},
    Xrot::AbstractVecOrMat{T}, eigval::Vector{T}, logdetV2::T)
    new(Yrot, Xrot, eigval, logdetV2)
  end
end

"""
Default constructor of a `TwoVarCompVariateRotate` instance.
"""
function TwoVarCompVariateRotate(
  Yrot::AbstractVecOrMat,
  Xrot::AbstractVecOrMat,
  eigval::Vector,
  logdetV2::Real)

  TwoVarCompVariateRotate{eltype(Yrot), typeof(Yrot), typeof(Xrot)}(Yrot, Xrot,
    eigval, logdetV2)
end

"""
Constructor of a `TwoVarCompVariateRotate` instance from a
`VarianceComponentVariate` instance.
"""
function TwoVarCompVariateRotate{T <: AbstractFloat}(
  vcobs::VarianceComponentVariate{T, 2}
  )

  # (generalized)-eigendecomposition of (V1, V2)
  if isa(vcobs.V[2], UniformScaling) ||
    (isdiag(vcobs.V[2]) && vecnorm(diag(vcobs.V[2]) - one(T)) < 1.0e-8)
    F = eigfact(Symmetric(vcobs.V[1]))
    deval = convert(Vector{T}, F.values)
    U = convert(Matrix{T}, F.vectors)
    logdetV2 = zero(T)
  else
    F = eigfact(Symmetric(vcobs.V[1]), Symmetric(vcobs.V[2]))
    deval = convert(Vector{T}, F.values)
    U = convert(Matrix{T}, F.vectors)
    logdetV2 = convert(T, logdet(vcobs.V[2]))
  end
  # corect negative eigenvalues due to roundoff error
  map!(x -> max(x, zero(T)), deval)::Vector{T}
  # rotate responses
  Yrot = At_mul_B(U, vcobs.Y)
  Xrot = isempty(vcobs.X) ? Array{T}(size(Yrot, 1), 0) : At_mul_B(U, vcobs.X)
  # output
  TwoVarCompVariateRotate(Yrot, Xrot, deval, logdetV2)
end

"""
Construct a `VarianceComponentModel` instance from a `VarianceComponentVariate`
instance. `B` is initialized to zero; `Σ` is initialized to a tupe of identity
matrices.
"""
function VarianceComponentModel{T, M}(vcobs::VarianceComponentVariate{T, M})
  p, d, m = size(vcobs.X, 2), size(vcobs.Y, 2), length(vcobs.V)
  B = zeros(T, p, d)
  Σ = ntuple(x -> eye(T, d), m)::NTuple{M, Matrix{T}}
  VarianceComponentModel(B, Σ)
end

"""
Construct a `VarianceComponentModel` instance from a `TwoVarCompVariateRotate`
instance.
"""
function VarianceComponentModel(vcobsrot::TwoVarCompVariateRotate)
  p, d = size(vcobsrot.Xrot, 2), size(vcobsrot.Yrot, 2)
  T = eltype(vcobsrot)
  B = zeros(T, p, d)
  Σ = (eye(T, d), eye(T, d))
  VarianceComponentModel(B, Σ)
end

Base.eltype(vcm::VarianceComponentModel) = Base.eltype(vcm.B)
Base.eltype(vcobs::VarianceComponentVariate) = Base.eltype(vcobs.Y)
Base.eltype(vcmrot::TwoVarCompModelRotate) = Base.eltype(vcmrot.Brot)
Base.eltype(vcobsrot::TwoVarCompVariateRotate) = Base.eltype(vcobsrot.Yrot)
# Dimension of response, d
length(vcm::VarianceComponentModel) = size(vcm.Σ[1], 1)
length(vcobs::VarianceComponentVariate) = size(vcobs.Y, 2)
length(vcmrot::TwoVarCompModelRotate) = length(vcmrot.eigval)
length(vcobsrot::TwoVarCompVariateRotate) = size(vcobsrot.Yrot, 2)
# Size of response, (n, d)
size(vcobs::VarianceComponentVariate) = size(vcobs.Y)
size(vcobsrot::TwoVarCompVariateRotate) = size(vcobsrot.Yrot)
# Number of variance components, m
nvarcomps(vcm::VarianceComponentModel) = length(vcm.Σ)
nvarcomps(vcobs::VarianceComponentVariate) = length(vcobs.V)
nvarcomps(vcmrot::TwoVarCompModelRotate) = 2
nvarcomps(vcobsrot::TwoVarCompVariateRotate) = 2
# Number of mean parameters, p * d
nmeanparams(vcm::VarianceComponentModel) = length(vcm.B)
nmeanparams(vcmrot::TwoVarCompModelRotate) = length(vcmrot.Brot)
nmeanparams(vcobs::VarianceComponentVariate) = size(vcobs.X, 2) * size(vcobs.Y, 2)
nmeanparams(vcobsrot::TwoVarCompVariateRotate) =
  size(vcobsrot.Xrot, 2) * size(vcobsrot.Yrot, 2)
# Number of free parameters in Cholesky factors, m * d * (d + 1) / 2
nvarparams(vcm::Union{VarianceComponentModel, TwoVarCompModelRotate}) =
  nvarcomps(vcm) * binomial(length(vcm) + 1, 2)
# Number of model parameters, p * d + m * d * (d + 1) / 2
nparams(vcm::Union{VarianceComponentModel, TwoVarCompModelRotate}) =
  nmeanparams(vcm) + nvarparams(vcm)

"""
    cov!(C, vcm, vcobs)

Calculate the `nd x nd` covariance matrix of a variance component model at an
observation and over-write `C`.
"""
function cov!(
  C::AbstractMatrix,
  vcm::VarianceComponentModel,
  vcobs::VarianceComponentVariate
  )

  fill!(C, zero(eltype(vcm)))
  for i in 1:length(vcm.Σ)
    kronaxpy!(vcm.Σ[i], vcobs.V[i], C)
  end
  return C
end

"""
    cov(vcm, vcobs)

Calculate the `nd x nd` covariance matrix of a variance component model at an
observation.
"""
function cov(
  vcm::VarianceComponentModel,
  vcobs::VarianceComponentVariate
  )

  n, d = size(vcobs)
  C = zeros(eltype(vcm), n * d, n * d)
  cov!(C, vcm, vcobs)
end

"""
    mean!(μ, vcm, vcobs)

Calculate the `n x d` mean matrix of a variance component model at an observation
and over-write `μ`.
"""
function mean!(
  μ::AbstractMatrix,
  vcm::VarianceComponentModel,
  vcobs::VarianceComponentVariate
  )

  A_mul_B!(μ, vcobs.X, vcm.B)
end

"""
    mean(vcm)

Calculate the `n x d` mean matrix of a variance component model at an observation.
"""
function mean(vcm::VarianceComponentModel, vcobs::VarianceComponentVariate)
  μ = zeros(eltype(vcobs), size(vcobs.Y))
  mean!(μ, vcm, vcobs)
end

function residual(vcm::VarianceComponentModel, vcobs::VarianceComponentVariate)
  return isempty(vcobs.X)? vcobs.Y : vcobs.Y - vcobs.X * vcm.B
end

function residual(vcmrot::TwoVarCompModelRotate, vcobsrot::TwoVarCompVariateRotate)
  if isempty(vcmrot.Brot)
    return vcobsrot.Yrot * vcmrot.eigvec
  else
    return vcobsrot.Yrot * vcmrot.eigvec - vcobsrot.Xrot * vcmrot.Brot
  end
end

# utilities for multivariate calculus
include("multivariate_calculus.jl")
# source for fitting models with 2 variance components
include("two_variance_component.jl")

end # VarianceComponentModels

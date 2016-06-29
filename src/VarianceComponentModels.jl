module VarianceComponentModels

using MathProgBase, Ipopt, KNITRO#, Mosek
import Base: eltype, length, size
export VarianceComponentModel, VarianceComponentVariate,
  TwoVarCompModelRotate, TwoVarCompVariateRotate, residual,
  nmeanparams, nvarparams, nparams

# M is number of variance components.
type VarianceComponentModel{T <: AbstractFloat, M,
  BT <: AbstractVecOrMat, ΣT <: AbstractMatrix}
  # model parameters
  B::BT
  Σ::NTuple{M, ΣT}
  function VarianceComponentModel(B::AbstractVecOrMat{T},
    Σ::NTuple{M, AbstractMatrix{T}})
    new(B, Σ)
  end
end

"""
Constructor for `VarianceComponentModel` instance.
"""
function VarianceComponentModel{M}(
  B::AbstractVecOrMat,
  Σ::NTuple{M, AbstractMatrix}
  )

  VarianceComponentModel{eltype(B), M, typeof(B), eltype(Σ)}(B, Σ)
end

function VarianceComponentModel{M}(Σ::NTuple{M, AbstractMatrix})
  B = zeros(eltype(Σ[1]), 0, size(Σ[1], 1))
  VarianceComponentModel{eltype(B), M, typeof(B), eltype(Σ)}(B, Σ)
end

"""
`eigval, eigvec = eig(Σ[1], Σ[2])` and `Brot = B * eigvec`.
"""
type TwoVarCompModelRotate{T <: AbstractFloat, BT <: AbstractVecOrMat}
  Brot::BT
  eigval::Vector{T}
  eigvec::Matrix{T}
  logdetΣ2::T
  function TwoVarCompModelRotate(Brot::AbstractVecOrMat{T},
    eigval::Vector{T}, eigvec::Matrix{T}, logdetΣ2::T)
    new(Brot, eigval, eigvec, logdetΣ2)
  end
end

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
function TwoVarCompModelRotate{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2}
  )

  # generalized eigenvalue decomposition of (Σ1, Σ2)
  λ, Φ = eig(vcm.Σ[1], vcm.Σ[2])
  # correct negative eigenvalues due to roundoff
  map!(x -> max(x, zero(T)), λ)
  Brot = isempty(vcm.B) ? Array{T}(0, size(vcm.Σ[1], 1)) : vcm.B * Φ
  logdetΣ2 = logdet(vcm.Σ[2])
  TwoVarCompModelRotate(Brot, λ, Φ, logdetΣ2)
end

immutable VarianceComponentVariate{T <: AbstractFloat, M,
  YT <: AbstractVecOrMat, XT <: AbstractVecOrMat, VT <: AbstractMatrix}
  # data
  Y::YT
  X::XT
  V::NTuple{M, VT}
  function VarianceComponentVariate(Y::AbstractVecOrMat{T},
    X::AbstractVecOrMat{T}, V::NTuple{M, AbstractMatrix{T}})
    new(Y, X, V)
  end
end

function VarianceComponentVariate{M}(
  Y::AbstractVecOrMat,
  X::AbstractVecOrMat,
  V::NTuple{M, AbstractMatrix}
  )

  VarianceComponentVariate{eltype(Y), M, typeof(Y), typeof(X), eltype(V)}(Y, X, V)
end

function VarianceComponentVariate{M}(
  Y::AbstractVecOrMat,
  V::NTuple{M, AbstractMatrix}
  )

  X = zeros(eltype(Y), size(Y, 1), 0)
  VarianceComponentVariate{eltype(Y), M, typeof(Y), typeof(X), eltype(V)}(Y, X, V)
end

"""
`eigval, eigvec = eig(Σ[1], Σ[2])`, `Yrot = eigvec * Y`, and 'Xrot = eigvec * X'.
"""
immutable TwoVarCompVariateRotate{T <: AbstractFloat, YT <: AbstractVecOrMat,
  XT <: AbstractVecOrMat}
  Yrot::YT
  Xrot::XT
  eigval::Vector{T}
  logdetV2::T
  function TwoVarCompVariateRotate(Yrot::AbstractVecOrMat{T},
    Xrot::AbstractVecOrMat{T}, eigval::Vector{T}, logdetV2::T)
    new(Yrot, Xrot, eigval, logdetV2)
  end
end

function TwoVarCompVariateRotate(
  Yrot::AbstractVecOrMat,
  Xrot::AbstractVecOrMat,
  eigval::Vector,
  logdetV2)

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
    deval, U = eig(vcobs.V[1])
    logdetV2 = zero(T)
  else
    deval, U = eig(vcobs.V[1], vcobs.V[2])
    logdetV2 = logdet(vcobs.V[2])
  end
  # corect negative eigenvalues due to roundoff error
  map!(x -> max(x, zero(eltype(vcobs))), deval)
  # rotate responses
  Yrot = At_mul_B(U, vcobs.Y)
  Xrot = isempty(vcobs.X) ? Array{T}(size(Yrot, 1), 0) : At_mul_B(U, vcobs.X)
  # output
  TwoVarCompVariateRotate(Yrot, Xrot, deval, logdetV2)
end

"""
Construct a `VarianceComponentModel` instance from a `VarianceComponentVariate`
instance.
"""
function VarianceComponentModel(vcobs::VarianceComponentVariate)
  p, d, m = size(vcobs.X, 2), size(vcobs.Y, 2), length(vcobs.V)
  B = zeros(eltype(vcobs), p, d)
  #Σ = (ntuple(x -> zeros(eltype(vcobs), d, d), m - 1)..., eye(eltype(vcobs), d))
  Σ = ntuple(x -> eye(eltype(vcobs), d), m)
  VarianceComponentModel{eltype(B), m, typeof(B), eltype(Σ)}(B, Σ)
end

Base.eltype(vcm::VarianceComponentModel) = Base.eltype(vcm.B)
Base.eltype(vcobs::VarianceComponentVariate) = Base.eltype(vcobs.Y)
Base.eltype(vcmrot::TwoVarCompModelRotate) = Base.eltype(vcmrot.Brot)
Base.eltype(vcobsrot::TwoVarCompVariateRotate) = Base.eltype(vcobsrot.Yrot)
# Dimension of response, d
length(vcm::VarianceComponentModel) = size(vcm.Σ[1], 1)
length(vcobs::VarianceComponentVariate) = size(vcobs.Y, 2)
length(vcmrot::TwoVarCompModelRotate) = length(vcmrot.eigval)
length(vcobsrot::TwoVarCompVariateRotate) = length(vcobsrot.eigval)
# Size of response, (n, d)
size(vcobs::VarianceComponentVariate) = size(vcobs.Y)
size(vcobsrot::TwoVarCompVariateRotate) = size(vcobsrot.Yrot)
# Number of variance components, m
nvarcomps(vcm::VarianceComponentModel) = length(vcm.Σ)
nvarcomps(vcobs::VarianceComponentVariate) = length(vcm.V)
nvarcomps(vcmrot::TwoVarCompModelRotate) = 2
nvarcomps(vcobsrot::TwoVarCompVariateRotate) = 2
# Number of mean parameters, p * d
nmeanparams(vcm::VarianceComponentModel) = length(vcm.B)
nmeanparams(vcmrot::TwoVarCompModelRotate) = length(vcmrot.Brot)
nmeanparams(vcobs::VarianceComponentVariate) = size(vcobs.X, 2) * size(vcobs.Y, 2)
nmeanparams(vcobsrot::TwoVarCompVariateRotate) =
  size(vcobsrot.Xrot, 2) * size(vcobs.Yrot, 2)
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
    kronaxpy!(C, vcm.Σ[i], vcobs.V[i])
  end
end

"""
    cov(vcm)

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
and over-write μ.
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
  μ = zeros(eltype(vcm), size(vcm.B))
  mean!(μ, vcm, vcobs)
end

function residual(vcm::VarianceComponentModel, vcobs::VarianceComponentVariate)
  resid = isempty(vcobs.X)? vcobs.Y : vcobs.Y - vcobs.X * vcm.B
end

function residual(vcmrot::TwoVarCompModelRotate, vcobsrot::TwoVarCompVariateRotate)
  if isempty(vcmrot.Brot)
    resid = vcobsrot.Yrot * vcmrot.eigvec
  else
    resid = vcobsrot.Yrot * vcmrot.eigvec - vcobsrot.Xrot * vcmrot.Brot
  end
  resid
end

# utilities for multivariate calculus
include("multivariate_calculus.jl")
# source for fitting models with 2 variance components
include("two_variance_component.jl")

end # VarianceComponentModels

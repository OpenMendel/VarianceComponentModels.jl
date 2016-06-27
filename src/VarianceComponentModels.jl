module VarianceComponentModels

using MathProgBase, Ipopt, KNITRO#, Mosek
import Base: eltype, length, size
export VarianceComponentModel, VarianceComponentVariate,
  TwoVarCompModelRotate, TwoVarCompVariateRotate

# M is number of variance components.
type VarianceComponentModel{T <: AbstractFloat, M}
  # model parameters
  B::AbstractVecOrMat{T}
  Σ::NTuple{M, AbstractMatrix{T}}
end

"""
`eigval, eigvec = eig(Σ[1], Σ[2])` and `Brot = B * eigvec`.
"""
type TwoVarCompModelRotate{T}
  Brot::AbstractVecOrMat{T}
  eigval::AbstractVecOrMat{T}
  eigvec::AbstractMatrix{T}
  logdetΣ2::T
end

"""
Constructor of a `TwoVarCompModelRotate` instance from a
`VarianceComponentModel` instance.
"""
function TwoVarCompModelRotate{T}(vcm::VarianceComponentModel{T, 2})
  # generalized eigenvalue decomposition of (Σ1, Σ2)
  λ, Φ = eig(vcm.Σ[1], vcm.Σ[2])
  # correct negative eigenvalues due to roundoff
  map!(x -> max(x, zero(T)), λ)
  Brot = isempty(vcm.B) ? T[] : A_mul_B(vcm.B, Φ)
  logdetΣ2 = logdet(vcm.Σ[2])
  TwoVarCompModelRotate(Brot, λ, Φ, logdetΣ2)
end

immutable VarianceComponentVariate{T <: AbstractFloat, M}
  # data
  Y::AbstractVecOrMat{T}
  X::AbstractVecOrMat{T}
  V::NTuple{M, AbstractMatrix{T}}
end

"""
`eigval, eigvec = eig(Σ[1], Σ[2])`, `Yrot = eigvec * Y`, and 'Xrot = eigvec * X'.
"""
immutable TwoVarCompVariateRotate{T}
  Yrot::AbstractVecOrMat{T}
  Xrot::AbstractVecOrMat{T}
  eigval::AbstractVecOrMat{T}
  logdetV2::T
end

"""
Constructor of a `TwoVarCompVariateRotate` instance from a
`VarianceComponentVariate` instance.
"""
function TwoVarCompVariateRotate{T <: Real}(
  vcobs::VarianceComponentVariate{T, 2}
  )

  # (generalized)-eigendecomposition of (V1, V2)
  if isa(vcobs.V[2], UniformScaling) ||
    (isdiag(vcobs.V[2]) && vecnorm(diag(vcobs.V[2]) - 1) < 1.0e-8)
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
  Xrot = isempty(vcobs.X) ? T[] : At_mul_B(U, vcobs.X)
  # output
  TwoVarCompVariateRotate(Yrot, Xrot, deval, logdetV2)
end

Base.eltype(vcm::VarianceComponentModel) = Base.eltype(vcm.B)
Base.eltype(vcobs::VarianceComponentVariate) = Base.eltype(vcobs.Y)
Base.eltype(vcmrot::TwoVarCompModelRotate) = Base.eltype(vcmrot.Brot)
Base.eltype(vcobsrot::TwoVarCompVariateRotate) = Base.eltype(vcobsrot.Yrot)
# Dimension of response, d
length(vcm::VarianceComponentModel) = size(vcm.Σ[1], 1)
length(vcobs::VarianceComponentVariate) = size(vcobs.Y, 2)
length(vcmrot::TwoVarCompModelRotate) = size(vcmrot.eigval, 1)
length(vcobs::TwoVarCompVariateRotate) = size(vcobs.Yrot, 2)
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
  μ::AbstractVecOrMat,
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
  if isempty(vcobsrot.Xrot)
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

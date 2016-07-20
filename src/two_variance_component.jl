import Base.gradient

export heritability,
  logpdf,
  gradient!, gradient,
  fisher_Σ!, fisher_Σ, fisher_B!, fisher_B,
  mle_fs!, mle_mm!,
  fit_mle!, fit_reml!,
  update_meanparam!

#---------------------------------------------------------------------------#
# Evaluate log-pdf
#---------------------------------------------------------------------------#

"""
    logpdf(vcmrot, vcobsrot)

Calculate log-pdf of a [`TwoVarCompVariateRotate`](@ref) instance under a
[`TwoVarCompModelRotate`](@ref).
"""
function logpdf{
  T1 <: TwoVarCompModelRotate,
  T2 <: TwoVarCompVariateRotate,
  T3 <: VarianceComponentAuxData}(
  vcmrot::T1,
  vcobsrot::T2,
  vcaux::T3 = VarianceComponentAuxData(vcobsrot)
  )

  n, d = size(vcobsrot.Yrot, 1), size(vcobsrot.Yrot, 2)
  T = eltype(vcmrot)
  zeroT, oneT = zero(T), one(T)
  residual!(vcaux.res, vcmrot, vcobsrot)
  # evaluate 2(log-likehood)
  objval = convert(T, - n * d * log(2π) - d * vcobsrot.logdetV2 - n * vcmrot.logdetΣ2)
  tmp, λj = zeroT, zeroT
  @inbounds for j in 1:d
    λj = vcmrot.eigval[j]
    @simd for i in 1:n
      tmp = oneT / (vcobsrot.eigval[i] * λj + oneT)
      objval += log(tmp) - tmp * vcaux.res[i, j]^2
    end
  end
  objval /= 2
end

"""
    logpdf(vcm, vcobsrot)

Calculate log-pdf of a [`TwoVarCompVariateRotate`](@ref) instance under a
[`VarianceComponentModel`](@ref).
"""
function logpdf{
  T1 <: VarianceComponentModel,
  T2 <: Union{TwoVarCompVariateRotate, VarianceComponentVariate},
  T3 <: VarianceComponentAuxData}(
  vcm::T1,
  vcobs::T2,
  vcaux::T3 = VarianceComponentAuxData(vcobs)
  )

  logpdf(TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs), vcaux)
end

"""
Calculate log-pdf of an array of variance component instances.
"""
function logpdf{
  T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate},
  T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate},
  T3 <: VarianceComponentAuxData}(
  vcm::T1,
  vcdata::Array{T2},
  vcaux::Array{T3} = map(x -> VarianceComponentAuxData(x), vcdata)
  )

  mapreduce(x -> logpdf(vcm, x...), +, zip(vcdata, vcaux))
end

#---------------------------------------------------------------------------#
# Evaluate gradient
#---------------------------------------------------------------------------#

"""
    gradient!(∇, vcmrot, vcobsrot)

Evaluate gradient of `Σ` at `vcmrot.Σ` and overwrite `∇`.

# Input
- `∇::Vector`: gradient vector.
- `vcmrot`: *rotated* two variance component model [`TwoVarCompModelRotate`](@ref).
- `vcobsrot`: *rotated* two variance component data [`TwoVarCompVariateRotate`](@ref).

# Output
- `∇`: gradient vector at `Σ = (Σ[1], Σ[2])`.

# TODO
- optimize computation
"""
function gradient!{T <: AbstractFloat}(
  ∇::AbstractVector{T},
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T},
  vcaux::VarianceComponentAuxData = VarianceComponentAuxData(vcobsrot)
  )

  n, d = size(vcobsrot.Yrot, 1), size(vcobsrot.Yrot, 2)
  zeroT, oneT = zero(T), one(T)
  residual!(vcaux.res, vcmrot, vcobsrot)
  # evaluate gradient with respect to Σ[1], Σ[2]
  m1diag = zeros(T, d)
  m2diag = zeros(T, d)
  tmp, λj = zeroT, zeroT
  @inbounds for j in 1:d
    λj = vcmrot.eigval[j]
    @simd for i in 1:n
      tmp = oneT / (vcobsrot.eigval[i] * λj + oneT)
      vcaux.res[i, j] *= tmp
      m1diag[j] += vcobsrot.eigval[i] * tmp
      m2diag[j] += tmp
    end
  end
  N2 = At_mul_B(vcaux.res, vcaux.res)
  scale!(sqrt(vcobsrot.eigval), vcaux.res)
  N1 = At_mul_B(vcaux.res, vcaux.res)
  @inbounds for j in 1:d
    N1[j, j] -= m1diag[j]
    N2[j, j] -= m2diag[j]
  end
  N1 = vcmrot.eigvec * N1 * vcmrot.eigvec'
  N2 = vcmrot.eigvec * N2 * vcmrot.eigvec'
  ∇[1:d^2] = N1[:]
  ∇[(d^2 + 1):(2d^2)] = N2[:]
  scale!(∇, convert(T, 0.5))
end # function gradient!

"""

    gradient(vcmrot, vcobsrot)

Evaluate gradient of `Σ` at `vcmrot.Σ` and overwrite `∇`.

# Input
- `vcmrot`: *rotated* two variance component model [`TwoVarCompModelRotate`](@ref).
- `vcobsrot`: *rotated* two variance component data [`TwoVarCompVariateRotate`](@ref).

# Output
- `∇`: gradient vector at `Σ = (Σ[1], Σ[2])`.

# TODO
- optimize computation
"""
function gradient{T <: AbstractFloat}(
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T},
  vcaux::VarianceComponentAuxData = VarianceComponentAuxData(vcobsrot)
  )

  d = length(vcmrot)
  ∇ = zeros(T, 2d^2)
  gradient!(∇, vcmrot, vcobsrot, vcaux)
end

"""
    gradient!(∇, vcm, vcobsrot)

Evaluate gradient of `Σ` at `vcmrot.Σ` and overwrite `∇`.

# Input
- `∇::Vector`: gradient vector.
- `vcmrot`: two variance component model [`VarianceComponentModel`](@ref).
- `vcobsrot`: *rotated* two variance component data [`TwoVarCompVariateRotate`](@ref).

# Output
- `∇`: gradient vector at `Σ = (Σ[1], Σ[2])`.

# TODO
- optimize computation
"""
function gradient!{T <: AbstractFloat}(
  ∇::AbstractVector{T},
  vcm::VarianceComponentModel{T, 2},
  vcobs::Union{TwoVarCompVariateRotate{T}, VarianceComponentVariate{T, 2}},
  vcaux::VarianceComponentAuxData = VarianceComponentAuxData(vcobs)
  )

  gradient!(∇, TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs), vcaux)
end


"""
    gradient(vcm, vcobsrot)

Evaluate gradient of `Σ` at `vcmrot.Σ` and overwrite `∇`.

# Input
- `vcmrot`: two variance component model [`VarianceComponentModel`](@ref).
- `vcobsrot`: *rotated* two variance component data [`TwoVarCompVariateRotate`](@ref).

# Output
- `∇`: gradient vector at `Σ = (Σ[1], Σ[2])`.

# TODO
- optimize computation
"""
function gradient{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobs::Union{TwoVarCompVariateRotate{T}, VarianceComponentVariate{T, 2}},
  vcaux::VarianceComponentAuxData = VarianceComponentAuxData(vcobs)
  )

  gradient(TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs), vcaux)
end

function gradient!{
  T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate},
  T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate},
  T3 <: VarianceComponentAuxData}(
  ∇::AbstractVector,
  vcm::T1,
  vcdata::Array{T2},
  vcaux::Array{T3} = map(x -> VarianceComponentAuxData(x), vcdata)
  )

  copy!(∇, gradient(vcm, vcdata, vcaux))
end


function gradient{
  T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate},
  T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate},
  T3 <: VarianceComponentAuxData}(
  vcm::T1,
  vcdata::Array{T2},
  vcaux::Array{T3} = map(x -> VarianceComponentAuxData(x), vcdata)
  )

  mapreduce(x -> gradient(vcm, x...), +, zip(vcdata, vcaux))
end

#---------------------------------------------------------------------------#
# Evaluate Fisher information matrix with respect to Σ
#---------------------------------------------------------------------------#

"""
    fisher_Σ!(H, vcmrot, vcobsrot)

Calculate Fisher information matrix at `Σ = (Σ[1], Σ[2])` and overwrite `H`
based on *rotated* two variance component model `vcmrot` and *rotated* two
variance component data `vcobsrot`.

# Input
- `H`: Hessian matrix.
- `vcmrot::TwoVarCompModelRotate`: *rotated* two variance component model.
- `vcobsrot::TwoVarCompVariateRotate`: *rotated* two variance component data.

# Output
- `H`: Fisher information matrix at `Σ = (Σ[1], Σ[2])`.
"""
function fisher_Σ!{T <: AbstractFloat}(
  H::AbstractMatrix{T},
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  d = length(vcmrot)
  zeroT, oneT = zero(T), one(T)
  fill!(H, zeroT)
  # evaluate Hessian with respect to Σ[1], Σ[2]
  C = zeros(T, d, d)
  Φ2 = kron(vcmrot.eigvec, vcmrot.eigvec)
  # (1, 1) block
  C[:] = T[mapreduce(
    x -> x^2 / (vcmrot.eigval[i] * x + oneT) / (vcmrot.eigval[j] * x + oneT),
    +, vcobsrot.eigval) for i=1:d, j=1:d]
  A_mul_Bt!(sub(H, 1:d^2, 1:d^2), scale(Φ2, vec(C)), Φ2)
  # (2, 1) block
  C[:] = T[mapreduce(
    x -> x / (vcmrot.eigval[i] * x + oneT) / (vcmrot.eigval[j] * x + oneT), +,
    vcobsrot.eigval) for i=1:d, j=1:d]
  A_mul_Bt!(sub(H, (d^2+1):(2d^2), 1:d^2), scale(Φ2, vec(C)), Φ2)
  # d-by-d (2, 2) block
  C[:] = T[mapreduce(
    x -> oneT / (vcmrot.eigval[i] * x + oneT) / (vcmrot.eigval[j] * x + oneT), +,
    vcobsrot.eigval) for i=1:d, j=1:d]
  A_mul_Bt!(sub(H, (d^2+1):(2d^2), (d^2+1):(2d^2)), scale(Φ2, vec(C)), Φ2)
  # copy to upper triangular part
  LinAlg.copytri!(H, 'L')
  scale!(H, convert(T, 0.5))
  # output
  return H
end # function fisher!

function fisher_Σ{T <: AbstractFloat}(
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  d = length(vcmrot.eigval)
  H = zeros(T, 2d^2, 2d^2)
  fisher_Σ!(H, vcmrot, vcobsrot)
end

function fisher_Σ!{T <: AbstractFloat}(
  H::AbstractMatrix{T},
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::Union{TwoVarCompVariateRotate{T}, VarianceComponentVariate{T, 2}}
  )

  fisher_Σ!(H, TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobsrot))
end

function fisher_Σ{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::Union{TwoVarCompVariateRotate{T}, VarianceComponentVariate{T, 2}}
  )

  fisher_Σ(TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobsrot))
end

function fisher_Σ!{
  T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate},
  T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate}}(
  H::AbstractMatrix,
  vcm::T1,
  vcobs::Array{T2}
  )

  copy!(H, fisher_Σ(vcm, vcobs))
end

function fisher_Σ{
  T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate},
  T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate}}(
  vcm::T1,
  vcobs::Array{T2}
  )

  mapreduce(x -> fisher_Σ(vcm, x), +, vcobs)
end

#---------------------------------------------------------------------------#
# Evaluate Fisher information matrix with respect to B
#---------------------------------------------------------------------------#

"""

    fisher_B!(H, vcmrot, vcobsrot)

Calculate Fisher information matrix of `B` and overwrite `H`
based on two variance component model `vcmrot` and *rotated* two
variance component data `vcobsrot`.

# Input
- `H`: Hessian matrix.
- `vcmrot::VarianceComponentModel{T, 2}`: two variance component model.
- `vcobsrot::TwoVarCompVariateRotate`: *rotated* two variance component data.

# Output
- `H`: Fisher information matrix at `B`.
"""
function fisher_B!{T <: AbstractFloat}(
  H::AbstractMatrix{T},
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T},
  vcaux::VarianceComponentAuxData = VarianceComponentAuxData(vcobsrot)
  )

  zeroT, oneT = zero(T), one(T)
  # working X
  fill!(vcaux.Xwork, zeroT)
  kronaxpy!(vcmrot.eigvec', vcobsrot.Xrot, vcaux.Xwork)
  # working weights
  fill!(vcaux.obswt, zeroT)
  kronaxpy!(vcobsrot.eigval, vcmrot.eigval, vcaux.obswt)
  @inbounds @simd for i in eachindex(vcaux.obswt)
    vcaux.obswt[i] = oneT / √(vcaux.obswt[i] + oneT)
  end
  # weighted least squares
  scale!(vcaux.obswt, vcaux.Xwork)
  At_mul_B!(H, vcaux.Xwork, vcaux.Xwork)
end

function fisher_B{T <: AbstractFloat}(
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  H = zeros(T, nmeanparams(vcmrot), nmeanparams(vcmrot))
  fisher_B!(H, vcmrot, vcobsrot)
end


function fisher_B!{T <: AbstractFloat}(
  H::AbstractMatrix{T},
  vcm::Union{VarianceComponentModel{T, 2}, TwoVarCompVariateRotate{T}},
  vcobs::Union{TwoVarCompVariateRotate{T}, VarianceComponentVariate{T, 2}},
  vcaux::VarianceComponentAuxData = VarianceComponentAuxData(vcobs)
  )

  fisher_B!(H, TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs), vcaux)
end

function fisher_B{
  T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate},
  T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate},
  T3 <: VarianceComponentAuxData}(
  vcm::T1,
  vcobs::T2,
  vcaux::T3 = VarianceComponentAuxData(vcobs)
  )

  H = zeros(eltype(vcm), nmeanparams(vcm), nmeanparams(vcm))
  fisher_B!(H, TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs), vcaux)
end

function fisher_B!{
  T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate},
  T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate},
  T3 <: VarianceComponentAuxData}(
  H::AbstractMatrix,
  vcm::T1,
  vcdata::Array{T2},
  vcaux::Array{T3} = map(x -> VarianceComponentAuxData(x), vcdata)
  )

  copy!(H, fisher_B(vcm, vcdata, vcaux))
end

function fisher_B{
  T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate},
  T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate},
  T3 <: VarianceComponentAuxData}(
  vcm::T1,
  vcdata::Array{T2},
  vcaux::Array{T3} = map(x -> VarianceComponentAuxData(x), vcdata)
  )

  mapreduce(x -> fisher_B(vcm, x...), +, zip(vcdata, vcaux))
end

#---------------------------------------------------------------------------#
# Update mean parameters by quadratic programming
#---------------------------------------------------------------------------#

"""
Compute the quadratic and linear parts of weighted least squares criterion.
"""
function suffstats_for_B{T <: AbstractFloat}(
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T},
  vcaux::VarianceComponentAuxData = VarianceComponentAuxData(vcobsrot)
  )

  zeroT, oneT = zero(T), one(T)
  # working weights
  fill!(vcaux.obswt, zeroT)
  kronaxpy!(vcmrot.eigval, vcobsrot.eigval, vcaux.obswt)
  @inbounds @simd for i in eachindex(vcaux.obswt)
    vcaux.obswt[i] = oneT / √(vcaux.obswt[i] + oneT)
  end
  # working X
  fill!(vcaux.Xwork, zeroT)
  kronaxpy!(vcmrot.eigvec', vcobsrot.Xrot, vcaux.Xwork)
  scale!(vcaux.obswt, vcaux.Xwork)
  # working y
  A_mul_B!(vcaux.res, vcobsrot.Yrot, vcmrot.eigvec)
  @inbounds @simd for i in eachindex(vcaux.ywork)
    vcaux.ywork[i] = vcaux.res[i] * vcaux.obswt[i]
  end
  # output
  return At_mul_B(vcaux.Xwork, vcaux.Xwork), At_mul_B(vcaux.Xwork, vcaux.ywork)
end

function suffstats_for_B{
  T1 <: TwoVarCompModelRotate,
  T2 <: TwoVarCompVariateRotate,
  T3 <: VarianceComponentAuxData}(
  vcmrot::T1,
  vcdatarot::Array{T2},
  vcaux::Array{T3}
  )

  mapreduce(x -> suffstats_for_B(vcmrot, x...), +,
    zip(vcdatarot, vcaux))
end

function +(
  x::Tuple{AbstractMatrix, AbstractVector},
  y::Tuple{AbstractMatrix, AbstractVector}
  )

  x[1] + y[1], x[2] + y[2]
end

"""
Update mean parameters `vcm.B` by quadratic programming.

## Input
- `vcm::VarianceComponentModel{T, 2}`: two variance component model
- `vcdatarot::TwoVarCompVariateRotate`: *rotated* two variance component data
- `qpsolver`: QP solver
- `Q`: quadratic part in the QP
- `c`: linear part in the QP

## Output
- `vcm.B`: updated mean parameters
"""
function update_meanparam!{
  T1 <: VarianceComponentModel,
  T2 <: TwoVarCompVariateRotate,
  T3 <: VarianceComponentAuxData}(
  vcm::T1,
  vcdatarot::Union{T2, Array{T2}},
  qpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver,
  vcaux::Union{T3, Array{T3}}
  )

  # quick return if there is no mean parameters
  isempty(vcm.B) && return vcm.B
  # rotate the model
  vcmrot = TwoVarCompModelRotate(vcm)
  # accumlate the quadratic and linear parts of QP
  Q, c = suffstats_for_B(vcmrot, vcdatarot, vcaux)
  # quadratic programming
  qpsol = quadprog(-c, Q, vcm.A, vcm.sense, vcm.b, vcm.lb, vcm.ub, qpsolver)
  if qpsol.status ≠ :Optimal
    println("Error in quadratic programming $(qpsol.status)")
  end
  copy!(vcm.B, qpsol.sol)
  return vcm.B
end

#---------------------------------------------------------------------------#
# Fisher scoring algorithm
#---------------------------------------------------------------------------#

"""
`TwoVarCompOptProb` holds data and model for solving two variance component
problem

# Fields
- `vcmodel`: [`VarianceComponentModel`](@ref)
- `qcdatarot`: [`TwoVarCompVariateRotate`](@ref)
- `qpsolver`: `MathProgBase.SolverInterface.AbstractMathProgSolver`
- `L`: working variable, a tuple of two Cholesky factors `(L[1], L[2])`
- `∇Σ`: working variable, gradient with respect to `(Σ[1], Σ[2])`
- `HΣ`: working variable, Hessian with respect to `(Σ[1], Σ[2])`
- `HL`: working variable, Hessian with respect to `(L[1], L[2])`
- `Q`: working variable, quadratic part of QP for updating `B`
- `c`: working variable, linear part of QP for updating `B`
"""
type TwoVarCompOptProb{
  T1 <: VarianceComponentModel,
  T2 <: TwoVarCompVariateRotate,
  T3 <: AbstractMatrix,
  T4 <: AbstractVector,
  T5 <: Union{VarianceComponentAuxData, AbstractArray}} <: MathProgBase.AbstractNLPEvaluator
  # variance component model and data
  vcmodel::T1
  vcdatarot::Union{T2, Array{T2}}
  # QP solver for updating B given Σ
  qpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver
  # intermediate variables
  L::NTuple{2, T3} # Cholesky factors
  ∇Σ::T4           # graident wrt (Σ1, Σ2)
  HΣ::T3           # Hessian wrt (Σ1, Σ2)
  HL::T3           # Hessian wrt (L1, L2)
  vcaux::T5
end

"""
    TwoVarCompOptProb(vcm, vcobsrot)
    TwoVarCompOptProb(vcm, vcobsrot, qpsolver)

Constructor of [`TwoVarCompOptProb`](@ref) from [`VarianceComponentModel`](@ref) `vcm`,
[`TwoVarCompVariateRotate`](@ref) `vcobsrot`, and input quadratic programming
sovler `qpsolver`.
"""
function TwoVarCompOptProb{
  T1 <: VarianceComponentModel,
  T2 <: TwoVarCompVariateRotate}(
  vcm::T1,
  vcdatarot::Union{T2, Array{T2}},
  qpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver
    = MathProgBase.defaultQPsolver
  )

  T = eltype(vcm)
  d, pd = length(vcm), nmeanparams(vcm)
  # number of optimization parameters in variance
  nvar = nvarparams(vcm)
  # allocate intermediate variables
  L  = (zeros(T, d, d), zeros(T, d, d))
  ∇Σ = zeros(T, 2d^2)       # graident wrt (Σ1, Σ2)
  HΣ = zeros(T, 2d^2, 2d^2) # Hessian wrt (Σ1, Σ2)
  HL = zeros(T, nvar, nvar) # Hessian wrt Ls
  if typeof(vcdatarot) <: AbstractArray
    vcaux = map(x -> VarianceComponentAuxData(x), vcdatarot)
    return TwoVarCompOptProb{typeof(vcm), eltype(vcdatarot), typeof(HΣ),
      typeof(∇Σ), typeof(vcaux)}(vcm, vcdatarot, qpsolver, L, ∇Σ, HΣ, HL, vcaux)
  else
    vcaux = VarianceComponentAuxData(vcdatarot)
    return TwoVarCompOptProb{typeof(vcm), typeof(vcdatarot), typeof(HΣ),
      typeof(∇Σ), typeof(vcaux)}(vcm, vcdatarot, qpsolver, L, ∇Σ, HΣ, HL, vcaux)
  end
end

"""
    optimparam_to_vcparam(dd, x)

Translate optimization variables `x` to variance component parmeters
`dd.vcmodel.Σ` in [`TwoVarCompOptProb`](@ref).
"""
function optimparam_to_vcparam!{T}(dd::TwoVarCompOptProb, x::Vector{T})

  d        = size(dd.L[1], 1)
  nvar     = nvarparams(dd.vcmodel)
  nvarhalf = div(nvar, 2)
  # variance parameter
  dd.L[1][trilind(d)] = x[1:nvarhalf]
  dd.L[2][trilind(d)] = x[(nvarhalf+1):end]
  # exponentiate diagonal entries
  @inbounds @simd for i in 1:d
    dd.L[1][i, i] = exp(dd.L[1][i, i])
    dd.L[2][i, i] = exp(dd.L[2][i, i])
  end
  # Σi = Li Li'
  A_mul_Bt!(dd.vcmodel.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(dd.vcmodel.Σ[2], dd.L[2], dd.L[2])
  # make sure the last variance component is pos. def.
  ϵ = convert(T, 1e-8)
  clamp_diagonal!(dd.vcmodel.Σ[2], ϵ, T(Inf))
end

function MathProgBase.initialize(dd::TwoVarCompOptProb,
  requested_features::Vector{Symbol})

  for feat in requested_features
    if !(feat in [:Grad, :Jac, :Hess])
      error("Unsupported feature $feat")
    end
  end
end # function MathProgBase.initialize

MathProgBase.features_available(dd::TwoVarCompOptProb) = [:Grad, :Jac, :Hess]
MathProgBase.eval_g(dd::TwoVarCompOptProb, g, x) = nothing
MathProgBase.jac_structure(dd::TwoVarCompOptProb) = Int[], Int[]
MathProgBase.eval_jac_g(dd::TwoVarCompOptProb, J, x) = nothing

function MathProgBase.eval_f{T}(dd::TwoVarCompOptProb, x::Vector{T})

  # update variance parameter from optim variable x
  optimparam_to_vcparam!(dd, x)
  # update mean parameters
  update_meanparam!(dd.vcmodel, dd.vcdatarot, dd.qpsolver, dd.vcaux)
  # evaluate profile log-pdf
  logpdf(dd.vcmodel, dd.vcdatarot, dd.vcaux)
end # function MathProgBase.eval_f

function MathProgBase.eval_grad_f{T}(
  dd::TwoVarCompOptProb,
  grad_f::Vector{T},
  x::Vector{T}
  )

  d        = size(dd.L[1], 1)
  nvar     = nvarparams(dd.vcmodel)
  nvarhalf = div(nvar, 2)
  # update variance parameter from optim variable x
  optimparam_to_vcparam!(dd, x)
  # update mean parameters
  update_meanparam!(dd.vcmodel, dd.vcdatarot, dd.qpsolver, dd.vcaux)
  # gradient wrt (Σ[1], Σ[2])
  gradient!(dd.∇Σ, dd.vcmodel, dd.vcdatarot, dd.vcaux)
  # chain rule for gradient wrt Cholesky factor
  chol_gradient!(sub(grad_f, 1:nvarhalf), dd.∇Σ[1:d^2], dd.L[1])
  chol_gradient!(sub(grad_f, (nvarhalf+1):nvar), dd.∇Σ[(d^2+1):end], dd.L[2])
  # chain rule for exponential of diagonal entries
  @inbounds for j in 1:d
    # linear index of diagonal entries of L1
    idx = 1 + (j - 1) * d - div((j - 1) * (j - 2), 2)
    grad_f[idx] *= dd.L[1][j, j]
    # linear index of diagonal entries of L2
    idx += nvarhalf
    grad_f[idx] *= dd.L[2][j, j]
  end
  return grad_f
end # function MathProgBase.eval_grad_f

function MathProgBase.hesslag_structure(dd::TwoVarCompOptProb)
  nvar = nvarparams(dd.vcmodel)
  # linear indices for variance parameters
  ind2sub((nvar, nvar), trilind(nvar))
end # function MathProgBase.hesslag_structure

function MathProgBase.eval_hesslag{T}(dd::TwoVarCompOptProb, H::Vector{T},
  x::Vector{T}, σ::T, μ::Vector{T})

  d        = size(dd.L[1], 1)
  nvar     = nvarparams(dd.vcmodel)
  nvarhalf = div(nvar, 2)
  # update variance parameter from optim variable x
  optimparam_to_vcparam!(dd, x)
  # update mean parameters
  update_meanparam!(dd.vcmodel, dd.vcdatarot, dd.qpsolver, dd.vcaux)
  # Hessian wrt (Σ1, Σ2)
  fisher_Σ!(dd.HΣ, dd.vcmodel, dd.vcdatarot)
  # chain rule for Hessian wrt Cholesky factor
  # only the lower left triangle
  # (1, 1) block
  chol_gradient!(sub(dd.HL, 1:nvarhalf, 1:nvarhalf),
    chol_gradient(dd.HΣ[1:d^2, 1:d^2], dd.L[1])',
    dd.L[1])
  # (2, 1) block
  chol_gradient!(sub(dd.HL, nvarhalf+1:nvar, 1:nvarhalf),
    chol_gradient(dd.HΣ[(d^2+1):(2d^2), 1:d^2], dd.L[1])',
    dd.L[2])
  # (2, 2) block
  chol_gradient!(sub(dd.HL, nvarhalf+1:nvar, nvarhalf+1:nvar),
    chol_gradient(dd.HΣ[(d^2+1):(2d^2), (d^2+1):(2d^2)], dd.L[2])',
    dd.L[2])
  # chain rule for exponential of diagonal entries
  @inbounds for j in 1:d
    # linear index of diagonal entries of L1
    idx = 1 + (j - 1) * d - div((j - 1) * (j - 2), 2)
    dd.HL[:, idx] *= dd.L[1][j, j]
    dd.HL[idx, :] *= dd.L[1][j, j]
    # linear index of diagonal entries of L2
    idx += nvarhalf
    dd.HL[:, idx] *= dd.L[2][j, j]
    dd.HL[idx, :] *= dd.L[2][j, j]
  end
  # output
  copy!(H, vech(dd.HL))
  scale!(H, -σ)
end

"""
    mle_fs!(vcmodel, vcdatarot; maxiter, solver, qpsolver, verbose)

Find MLE by Fisher scoring algorithm.

# Input
- `vcmodel`: two variane component model [`VarianceComponentModel`](@ref), with
 `vcmodel.B` and `vcmodel.Σ` used as starting point
- `vcdatarot`: rotated two varianec component data [`TwoVarCompVariateRotate`](@ref)

# Keyword
- `maxiter::Int`: maximum number of iterations, default is 1000
- `solver::Symbol`: backend nonlinear programming solver, `:Ipopt` (default) or `:Knitro`
- `qpsolver::Symbol`: backend quadratic programming solver, `:Ipopt` (default) or `:Gurobi` or `Mosek`
- `verbose::Bool`: display information

# Output
- `maxlogl`: log-likelihood at solution
- `vcmodel`: [`VarianceComponentModel`](@ref) with updated model parameters
- `Σse=(Σse[1],Σse[2])`: standard errors of estimate `Σ=(Σ[1],Σ[2])`
- `Σcov`: covariance matrix of estimate `Σ=(Σ[1],Σ[2])`
- `Bse`: standard errors of estimate `B`
- `Bcov`: covariance of estimate `B`
"""
function mle_fs!{T1<:VarianceComponentModel, T2<:TwoVarCompVariateRotate}(
  vcmodel::T1,
  vcdatarot::Union{T2, Array{T2}};
  maxiter::Integer = 1000,
  solver::Symbol = :Ipopt,
  qpsolver::Symbol = :Ipopt,
  verbose::Bool = true,
  )

  T = eltype(vcmodel)
  d = length(vcmodel)
  Ltrilind = trilind(d, d)
  # number of mean parameters
  nmean = nmeanparams(vcmodel)
  # number of optimization parameters in Cholesky factors
  nvar = nvarparams(vcmodel)
  nvarhalf = div(nvar, 2)
  # pre-allocate variables for optimization
  zeroT = convert(T, 0)
  # data for the optimization problem
  if qpsolver == :Ipopt
    qs = IpoptSolver(print_level = 0)
  elseif qpsolver == :Gurobi
    qs = GurobiSolver(OutputFlag = 0)
  elseif qpsolver == :Mosek
    qs = MosekSolver(MSK_IPAR_LOG = 0)
  end
  dd = TwoVarCompOptProb(vcmodel, vcdatarot, qs)

  # set up MathProgBase interface
  if solver == :Ipopt
    # see http://www.coin-or.org/Ipopt/documentation/documentation.html for IPOPT
    solver = IpoptSolver(
      hessian_approximation = "exact",
      tol = 1.0e-8, # default is 1.0e-8
      acceptable_tol = 1.0e-5, # default is 1.0e-6
      max_iter = maxiter, # default is 3000
      print_frequency_iter = 5, # default is 1
      print_level = verbose? 5 : 0,
      print_info_string = "yes",
      #nlp_scaling_method = "none",
      #derivative_test = "second-order",
      #linear_solver = "mumps",
      #linear_solver = "pardiso",
      mehrotra_algorithm = "yes",
      )
  elseif solver == :Mosek
    # see http://docs.mosek.com/7.0/capi/Parameters.html for Mosek options
    solver = MosekSolver(
      MSK_IPAR_INTPNT_MAX_ITERATIONS = maxiter,
      MSK_DPAR_INTPNT_NL_TOL_REL_GAP = 1.0e-8,
      MSK_IPAR_LOG = verbose? 10 : 0, # deafult value is 10
      #MSK_IPAR_OPTIMIZER = MSK_OPTIMIZER_NONCONVEX,
      #MSK_IPAR_LOG_NONCONVEX = 20,
      #MSK_IPAR_NONCONVEX_MAX_ITERATIONS = 100,
      #MSK_DPAR_INTPNT_NL_TOL_NEAR_REL = 1e8,
      #MSK_IPAR_LOG_CHECK_CONVEXITY = 1,
      #MSK_IPAR_INFEAS_PREFER_PRIMAL = MSK_OFF
      )
  elseif solver == :Knitro
    # see https://www.artelys.com/tools/knitro_doc/3_referenceManual/userOptions.html for Mosek options
    solver = KnitroSolver(
      KTR_PARAM_ALG = 0,
      KTR_PARAM_OUTLEV = verbose? 2 : 0,
      #KTR_PARAM_MAXCGIT = 5,
      #KTR_PARAM_SCALE = 0,
      #KTR_PARAM_HONORBNDS = 1,
      #KTR_PARAM_GRADOPT = 1,
      #KTR_PARAM_HESSOPT = 1,
      #KTR_PARAM_DERIVCHECK = 3,
      #KTR_PARAM_TUNER = 1,
      #KTR_PARAM_MAXTIMECPU = 5.0,
      #KTR_PARAM_BAR_MURULE = 4,
      #KTR_PARAM_BAR_FEASIBLE = 3,
      #KTR_PARAM_BAR_FEASMODETOL = 1.0e-8,
      #KTR_PARAM_BAR_SWITCHRULE = 3,
      #KTR_PARAM_BAR_PENCONS = 2,
      #KTR_PARAM_BAR_DIRECTINTERVAL = 0,
      #KTR_PARAM_BAR_MAXBACKTRACK = 10,
      #KTR_PARAM_BAR_MAXCROSSIT = 5,
      #KTR_PARAM_BAR_MAXREFACTOR = 5,
      #KTR_PARAM_BAR_REFINEMENT = 1,
      #KTR_PARAM_BAR_WATCHDOG = 1,
      )
  end
  m = MathProgBase.NonlinearModel(solver)
  # lower and upper bounds for variance parameter
  lb = zeros(nvar); fill!(lb, -T(Inf))
  ub = zeros(nvar); fill!(ub,  T(Inf))
  MathProgBase.loadproblem!(m, nvar, 0, lb, ub, T[], T[], :Max, dd)
  # start point
  copy!(dd.L[1], cholfact(vcmodel.Σ[1], :L, Val{true})[:L].data)
  copy!(dd.L[2], cholfact(vcmodel.Σ[2], :L, Val{true})[:L].data)
  # reparameterize diagonal entries to exponential
  @inbounds @simd for i in 1:d
    dd.L[1][i, i] = log(dd.L[1][i, i] + convert(T, 1e-8))
    dd.L[2][i, i] = log(dd.L[2][i, i] + convert(T, 1e-8))
  end
  x0 = [vech(dd.L[1]); vech(dd.L[2])]
  MathProgBase.setwarmstart!(m, x0)
  # optimize
  MathProgBase.optimize!(m)
  stat = MathProgBase.status(m)
  x = MathProgBase.getsolution(m)
  # retrieve variance component parameters
  optimparam_to_vcparam!(dd, x)
  # update mean parameters
  update_meanparam!(dd.vcmodel, dd.vcdatarot, dd.qpsolver, dd.vcaux)
  # update final objective value
  maxlogl = logpdf(vcmodel, vcdatarot, dd.vcaux)

  # standard errors
  Bcov = zeros(T, nmean, nmean)
  Bcov = inv(fisher_B(vcmodel, vcdatarot, dd.vcaux))
  Bse = similar(vcmodel.B)
  copy!(Bse, sqrt(diag(Bcov)))
  Σcov = inv(fisher_Σ(vcmodel, vcdatarot))
  Σse = (zeros(T, d, d), zeros(T, d, d))
  copy!(Σse[1], sqrt(diag(sub(Σcov, 1:d^2, 1:d^2))))
  copy!(Σse[2], sqrt(diag(sub(Σcov, d^2+1:2d^2, d^2+1:2d^2))))

  # output
  maxlogl, vcmodel, Σse, Σcov, Bse, Bcov
end # function mle_fs

#---------------------------------------------------------------------------#
# MM algorithm
#---------------------------------------------------------------------------#

function suffstats_for_Σ!(
  vcmrot::TwoVarCompModelRotate,
  vcobsrot::TwoVarCompVariateRotate,
  vcaux::VarianceComponentAuxData = VarianceComponentAuxData(vcobsrot)
  )

  T = eltype(vcmrot)
  n, d = size(vcobsrot.Yrot, 1), size(vcobsrot.Yrot, 2)
  A1, b1 = zeros(T, d, d), zeros(T, d)
  A2, b2 = zeros(T, d, d), zeros(T, d)
  zeroT, oneT = zero(T), one(T)
  λj, tmp     = zeroT, zeroT
  # (rotated) residual
  residual!(vcaux.res, vcmrot, vcobsrot)
  # sufficient statistics for Σ2
  @inbounds for j in 1:d
    λj = vcmrot.eigval[j]
    for i in 1:n
      tmp = oneT / (vcobsrot.eigval[i] * λj + oneT)
      vcaux.res[i, j] *= tmp
      b2[j] += tmp
    end
  end
  At_mul_B!(A2, vcaux.res, vcaux.res)
  # sufficient statistics for Σ1
  @inbounds for j in 1:d
    λj = vcmrot.eigval[j]
    for i in 1:n
      tmp = vcobsrot.eigval[i]
      vcaux.res[i, j] *= √tmp * λj
      b1[j] += tmp / (tmp * λj + oneT)
    end
  end
  At_mul_B!(A1, vcaux.res, vcaux.res)
  # output
  return A1, b1, A2, b2
end

function suffstats_for_Σ!{
  T1 <: TwoVarCompModelRotate,
  T2 <: TwoVarCompVariateRotate,
  T3 <: VarianceComponentAuxData}(
  vcmrot::T1,
  vcdatarot::Array{T2},
  vcaux::Array{T3} = map(x -> VarianceComponentAuxData(x), vcdatarot)
  )

  mapreduce(x -> suffstats_for_Σ!(vcmrot, x...), +, zip(vcdatarot, vcaux))
end

function +(
  x::Tuple{AbstractMatrix, AbstractVector, AbstractMatrix, AbstractVector},
  y::Tuple{AbstractMatrix, AbstractVector, AbstractMatrix, AbstractVector}
  )

  x[1] + y[1], x[2] + y[2], x[3] + y[3], x[4] + y[4]
end

function mm_update_Σ!{
  T1 <: VarianceComponentModel,
  T2 <: TwoVarCompVariateRotate,
  T3 <: Union{VarianceComponentAuxData, AbstractArray}}(
  vcm::T1,
  vcdatarot::Union{T2, Array{T2}},
  vcaux::T3
  )

  T, d = eltype(vcm), length(vcm)
  zeroT, oneT = zero(T), one(T)
  # eigen-decomposition of (vcm.Σ[1], vcm.Σ[2])
  vcmrot = TwoVarCompModelRotate(vcm)
  # sufficient statistics for updating Σ1, Σ2
  A1, b1, A2, b2 = suffstats_for_Σ!(vcmrot, vcdatarot, vcaux)
  @inbounds for j in 1:d
    b1[j] = √b1[j]
    b2[j] = √b2[j]
  end
  Φinv = inv(vcmrot.eigvec)
  # update Σ1
  scale!(b1, A1), scale!(A1, b1)
  storage = eigfact!(Symmetric(A1))
  @inbounds for i in 1:d
    storage.values[i] = storage.values[i] > zeroT ? √√storage.values[i] : zeroT
  end
  scale!(storage.vectors, storage.values)
  scale!(oneT ./ b1, storage.vectors)
  At_mul_B!(vcm.Σ[1], Φinv, storage.vectors)
  A_mul_Bt!(vcm.Σ[1], vcm.Σ[1], vcm.Σ[1])
  # update Σ2
  scale!(b2, A2), scale!(A2, b2)
  storage = eigfact!(Symmetric(A2))
  @inbounds for i in 1:d
    storage.values[i] = storage.values[i] > zeroT ? √√storage.values[i] : zeroT
  end
  scale!(storage.vectors, storage.values)
  scale!(oneT ./ b2, storage.vectors)
  At_mul_B!(vcm.Σ[2], Φinv, storage.vectors)
  A_mul_Bt!(vcm.Σ[2], vcm.Σ[2], vcm.Σ[2])
end

"""
    mle_mm!(vcmodel, vcdatarot; maxiter, qpsolver, verbose)

Find MLE by minorization-maximization (MM) algorithm.

# Input
- `vcmodel`: two variane component model [`VarianceComponentModel`](@ref), with
 `vcmodel.B` and `vcmodel.Σ` used as starting point
- `vcdatarot`: rotated two varianec component data [`TwoVarCompVariateRotate`](@ref)

# Keyword
- `maxiter::Int`: maximum number of iterations, default is 1000
- `qpsolver::Symbol`: backend quadratic programming solver, `:Ipopt` (default) or `:Gurobi` or `Mosek`
- `verbose::Bool`: display information

# Output
- `maxlogl`: log-likelihood at solution
- `vcmodel`: [`VarianceComponentModel`](@ref) with updated model parameters
- `Σse=(Σse[1],Σse[2])`: standard errors of estimate `Σ=(Σ[1],Σ[2])`
- `Σcov`: covariance matrix of estimate `Σ=(Σ[1],Σ[2])`
- `Bse`: standard errors of estimate `B`
- `Bcov`: covariance of estimate `B`

# Reference
- H. Zhou, L. Hu, J. Zhou, and K. Lange (2015)
  MM algorithms for variance components models.
  [http://arxiv.org/abs/1509.07426](http://arxiv.org/abs/1509.07426)
"""
function mle_mm!{T1 <: VarianceComponentModel, T2 <: TwoVarCompVariateRotate}(
  vcm::T1,
  vcdatarot::Union{T2, Array{T2}};
  maxiter::Integer = 10000,
  funtol::Real = convert(eltype(vcm), 1e-8),
  qpsolver::Symbol = :Ipopt,
  verbose::Bool = true,
  )

  T, d, pd = eltype(vcm), length(vcm), nmeanparams(vcm)
  # set up QP solver
  if qpsolver == :Ipopt
    qs = IpoptSolver(print_level = 0)
  elseif qpsolver == :Gurobi
    qs = GurobiSolver(OutputFlag = 0)
  elseif qpsolver == :Mosek
    qs = MosekSolver(MSK_IPAR_LOG = 0)
  end
  # allocate intermediate variables
  if typeof(vcdatarot) <: AbstractArray
    vcaux = map(x -> VarianceComponentAuxData(x), vcdatarot)
  else
    vcaux = VarianceComponentAuxData(vcdatarot)
  end
  # initial log-likelihood
  logl::T = logpdf(vcm, vcdatarot, vcaux)
  if verbose
    println()
    println("     MM Algorithm")
    println("  Iter      Objective  ")
    println("--------  -------------")
    @printf("%8.d  %13.e\n", 0, logl)
  end

  # MM loop
  for iter in 1:maxiter

    # update Σ
    mm_update_Σ!(vcm, vcdatarot, vcaux)
    # make sure the last variance component is pos. def.
    ϵ = convert(T, 1e-8)
    clamp_diagonal!(vcm.Σ[2], ϵ, T(Inf))

    # update mean parameters
    if pd > 0; update_meanparam!(vcm, vcdatarot, qs, vcaux); end

    # check convergence
    loglold = logl
    logl    = logpdf(vcm, vcdatarot, vcaux)
    if verbose
      if (iter <= 10) || (iter > 10 && iter % 10 == 0)
        @printf("%8.d  %13.e\n", iter, logl)
      end
    end
    if abs(logl - loglold) < funtol * (abs(logl) + one(T))
      break
    end
  end
  if verbose; println(); end

  # standard errors
  Bcov = inv(fisher_B(vcm, vcdatarot, vcaux))
  Bse = similar(vcm.B)
  copy!(Bse, sqrt(diag(Bcov)))
  Σcov = inv(fisher_Σ(vcm, vcdatarot))
  Σse = (zeros(T, d, d), zeros(T, d, d))
  copy!(Σse[1], sqrt(diag(sub(Σcov, 1:d^2, 1:d^2))))
  copy!(Σse[2], sqrt(diag(sub(Σcov, d^2+1:2d^2, d^2+1:2d^2))))

  # output
  logl, vcm, Σse, Σcov, Bse, Bcov
end # function mle_mm

#---------------------------------------------------------------------------#
# Estimation gateway
#---------------------------------------------------------------------------#

"""
    fit_mle!(vcmodel, vcdata; algo)

Find MLE of variane component model.

# Input
- `vcmodel`: two variane component model [`VarianceComponentModel`](@ref), with
 `vcmodel.B` and `vcmodel.Σ` used as starting point
- `vcdata`: two varianec component data [`VarianceComponentVariate`](@ref)

# Keyword
- `algo::Symbol`: algorithm, `:FS` (Fisher scoring) for `:MM`
(minorization-maximization algorithm)

# Output
- `maxlogl`: log-likelihood at solution
- `vcmodel`: [`VarianceComponentModel`](@ref) with updated model parameters
- `Σse=(Σse[1],Σse[2])`: standard errors of estimate `Σ=(Σ[1],Σ[2])`
- `Σcov`: covariance matrix of estimate `Σ=(Σ[1],Σ[2])`
- `Bse`: standard errors of estimate `B`
- `Bcov`: covariance of estimate `B`
"""
function fit_mle!{
  T1 <: VarianceComponentModel,
  T2 <: VarianceComponentVariate}(
  vcmodel::T1,
  vcdata::Union{T2, Array{T2}};
  algo::Symbol = :FS
  )

  # generalized-eigendecomposition and rotate data
  vcdatarot = TwoVarCompVariateRotate(vcdata)
  if algo == :FS
    return mle_fs!(vcmodel, vcdatarot)
  elseif algo == :MM
    return mle_mm!(vcmodel, vcdatarot)
  end
end

"""
    fit_reml!(vcmodel, vcdata; algo)

Find restricted MLE (REML) of variane component model.

# Input
- `vcmodel`: two variane component model [`VarianceComponentModel`](@ref), with
 `vcmodel.B` and `vcmodel.Σ` used as starting point
- `vcdata`: two varianec component data [`VarianceComponentVariate`](@ref)

# Keyword
- `algo::Symbol`: algorithm, `:FS` (Fisher scoring) for `:MM`
(minorization-maximization algorithm)

# Output
- `maxlogl`: log-likelihood at solution
- `vcmodel`: [`VarianceComponentModel`](@ref) with updated model parameters
- `Σse=(Σse[1],Σse[2])`: standard errors of estimate `Σ=(Σ[1],Σ[2])`
- `Σcov`: covariance matrix of estimate `Σ=(Σ[1],Σ[2])`
- `Bse`: standard errors of estimate `B`
- `Bcov`: covariance of estimate `B`
"""
function fit_reml!{
  T1 <: VarianceComponentModel,
  T2 <: VarianceComponentVariate}(
  vcmodel::T1,
  vcdata::Union{T2, Array{T2}};
  algo::Symbol = :FS
  )

  d = length(vcmodel)
  vcdatarot = TwoVarCompVariateRotate(vcdata)
  # residual assuming (Σ[1], Σ[2]) = (eye(d), eye(d))
  vcmtmp = deepcopy(vcmodel)
  copy!(vcmtmp.Σ[1], eye(d)), copy!(vcmtmp.Σ[2], eye(d))
  update_meanparam!(vcmtmp, vcdatarot)
  resrot = residual(TwoVarCompModelRotate(vcmtmp), vcdatarot)
  # use residuals as responses
  resdatarot = TwoVarCompVariateRotate(resrot, zeros(T, n, 0),
    vcdatarot.eigval, vcdatarot.logdetV2)
  resmodel = deepcopy(vcmodel)
  resmodel.B = zeros(T, 0, d)
  if algo == :FS
    _, _, Σse, Σcov, = mle_fs!(resmodel, resdatarot)
  elseif algo == :MM
    _, _, Σse, Σcov, = mle_mm!(resmodel, resdatarot)
  end

  # estimate mean parameters from covariance estimate
  copy!(vcmodel.Σ[1], resmodel.Σ[1]), copy!(vcmodel.Σ[2], resmodel.Σ[2])
  update_meanparam!(vcmodel, vcdatarot)

  # standard errors and covariance of mean parameters
  Bcov = inv(fisher_B(vcmodel, vcdatarot))
  Bse = similar(vcmodel.B)
  copy!(Bse, sqrt(diag(Bcov)))

  # output
  logpdf(vcmodel, vcdatarot), vcmodel, Σse, Σcov, Bse, Bcov
end

#---------------------------------------------------------------------------#
# Heritability estimation
#---------------------------------------------------------------------------#

"""
    heritability(Σ, Σcov, which=1)

Calcualte the heritability of each trait and its standard error from variance
component estimates and their covariance.

# Input
- `Σ=(Σ[1], ..., Σ[m])`: estimate of `m` `d x d` variance components.
- `Σcov`: `md^2 x md^2` covariance matrix of the variance component estiamte.
- `which`: indicator which is the additive genetic variance component.

# Output
- `h`: estimated heritability of `d` traits.
- `hse`: standard errors of the estiamted heritability.
"""
function heritability{T, N}(
  Σ::NTuple{N, AbstractMatrix{T}},
  Σcov::AbstractMatrix{T},
  which::Integer = 1
  )

  d = size(Σ[1], 1) # number of traits
  m = length(Σ)     # number of variance components
  @assert size(Σcov, 1) == m * d^2 "Dimensions don't match"
  h, hse = zeros(T, d), zeros(T, d)
  hgrad = zeros(T, m)
  for trait in 1:d
    σ2trait = [Σ[i][trait, trait] for i = 1:m]
    totalσ2trait = sum(σ2trait)
    # heritability
    h[trait] = σ2trait[which] / totalσ2trait
    # standard error of heritability by Delta method
    for i in 1:m
      if i == which
        hgrad[i] = (totalσ2trait - σ2trait[i]) / totalσ2trait^2
      else
        hgrad[i] = - σ2trait[i] / totalσ2trait^2
      end
    end
    dgidx = diagind(Σ[1])[trait]
    traitidx = dgidx:d^2:((m - 1) * d^2 + dgidx)
    hse[trait] = sqrt(dot(hgrad, sub(Σcov, traitidx, traitidx) * hgrad))
  end
  # output
  h, hse
end # function heritability

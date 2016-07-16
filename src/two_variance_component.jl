import Base.gradient

export heritability,
  logpdf, gradient!, gradient, fisher!, fisher, fisher_B!, fisher_B,
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
function logpdf{T <: AbstractFloat}(
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  n, d = size(vcobsrot.Yrot, 1), size(vcobsrot.Yrot, 2)
  nd = n * d
  zeroT, oneT = zero(T), one(T)
  res = residual(vcmrot, vcobsrot)
  # evaluate 2(log-likehood)
  objval = convert(T, - nd * log(2π) - d * vcobsrot.logdetV2 - n * vcmrot.logdetΣ2)
  tmp, λj = zeroT, zeroT
  @inbounds for j in 1:d
    λj = vcmrot.eigval[j]
    @simd for i in 1:n
      tmp = oneT / (vcobsrot.eigval[i] * λj + oneT)
      objval += log(tmp) - tmp * res[i, j]^2
    end
  end
  objval /= 2
end

"""
    logpdf(vcm, vcobsrot)

Calculate log-pdf of a [`TwoVarCompVariateRotate`](@ref) instance under a
[`VarianceComponentModel`](@ref).
"""
function logpdf{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  logpdf(TwoVarCompModelRotate(vcm), vcobsrot)
end

"""
    logpdf(vcm, vcobsrot)

Calculate log-pdf of a [`VarianceComponentVariate`](@ref) instance under a
[`VarianceComponentModel`](@ref).
"""
function logpdf{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::VarianceComponentVariate{T, 2}
  )

  logpdf(TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobsrot))
end

"""
Calculate log-pdf of an array of variance component instances.
"""
function logpdf{T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate}, T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate}}(
  vcm::T1,
  vcobs::Array{T2}
  )

  mapreduce(x -> logpdf(vcm, x), +, vcobs)
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
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  n, d = size(vcobsrot.Yrot, 1), size(vcobsrot.Yrot, 2)
  zeroT, oneT = zero(T), one(T)
  res = residual(vcmrot, vcobsrot)
  # evaluate gradient with respect to Σ[1], Σ[2]
  m1diag = zeros(T, d)
  m2diag = zeros(T, d)
  tmp, λj = zeroT, zeroT
  @inbounds for j in 1:d
    λj = vcmrot.eigval[j]
    @simd for i in 1:n
      tmp = oneT / (vcobsrot.eigval[i] * λj + oneT)
      res[i, j] *= tmp
      m1diag[j] += vcobsrot.eigval[i] * tmp
      m2diag[j] += tmp
    end
  end
  N2 = At_mul_B(res, res)
  scale!(sqrt(vcobsrot.eigval), res)
  N1 = At_mul_B(res, res)
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
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  d = length(vcmrot)
  ∇ = zeros(T, 2d^2)
  gradient!(∇, vcmrot, vcobsrot)
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
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  gradient!(∇, TwoVarCompModelRotate(vcm), vcobsrot)
end

"""
    gradient!(∇, vcm, vcobsrot)

Evaluate gradient of `Σ` at `vcmrot.Σ` and overwrite `∇`.

# Input
- `∇::Vector`: gradient vector.
- `vcm`: variance component model [`VarianceComponentModel`](@ref).
- `vcobs`: two variance component data [`VarianceComponentVariate`](@ref).

# Output
- `∇`: gradient vector at `Σ = (Σ[1], Σ[2])`.

# TODO
- optimize computation
"""
function gradient!{T <: AbstractFloat}(
  ∇::AbstractVector{T},
  vcm::VarianceComponentModel{T, 2},
  vcobs::VarianceComponentVariate{T, 2}
  )

  gradient!(∇, TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs))
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
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  gradient(TwoVarCompModelRotate(vcm), vcobsrot)
end

"""
    gradient(vcm, vcobsrot)

Evaluate gradient of `Σ` at `vcmrot.Σ` and overwrite `∇`.

# Input
- `vcm`: two variance component model [`VarianceComponentModel`](@ref).
- `vcobs`: two variance component data [`VarianceComponentModel`](@ref).

# Output
- `∇`: gradient vector at `Σ = (Σ[1], Σ[2])`.

# TODO
- optimize computation
"""
function gradient{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobs::VarianceComponentVariate{T, 2}
  )

  gradient(TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs))
end

function gradient!{T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate}, T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate}}(
  ∇::AbstractVector,
  vcm::T1,
  vcobs::Array{T2}
  )

  copy!(∇, gradient(vcm, vcobs))
end


function gradient{T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate}, T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate}}(
  vcm::T1,
  vcobs::Array{T2}
  )

  mapreduce(x -> gradient(vcm, x), +, vcobs)
end

#---------------------------------------------------------------------------#
# Evaluate Fisher information matrix with respect to Σ
#---------------------------------------------------------------------------#

"""
    fisher!(H, vcmrot, vcobsrot)

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
function fisher!{T <: AbstractFloat}(
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
end # function fisher!

function fisher{T <: AbstractFloat}(
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  d = length(vcmrot.eigval)
  H = zeros(T, 2d^2, 2d^2)
  fisher!(H, vcmrot, vcobsrot)
end

function fisher!{T <: AbstractFloat}(
  H::AbstractMatrix{T},
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  fisher!(H, TwoVarCompModelRotate(vcm), vcobsrot)
end

function fisher{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  fisher(TwoVarCompModelRotate(vcm), vcobsrot)
end

function fisher!{T <: AbstractFloat}(
  H::AbstractMatrix{T},
  vcm::VarianceComponentModel{T, 2},
  vcobs::VarianceComponentVariate{T, 2}
  )

  fisher!(H, TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs))
end

function fisher{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobs::VarianceComponentVariate{T, 2}
  )

  fisher(TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs))
end


function fisher!{T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate}, T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate}}(
  H::AbstractMatrix,
  vcm::T1,
  vcobs::Array{T2}
  )

  copy!(H, fisher(vcm, vcobs))
end

function fisher{T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate}, T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate}}(
  vcm::T1,
  vcobs::Array{T2}
  )

  mapreduce(x -> fisher(vcm, x), +, vcobs)
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
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  oneT = one(T)
  M = kron(vcmrot.eigvec, vcobsrot.Xrot')
  A_mul_Bt!(H, scale(M, oneT ./ (kron(vcobsrot.eigval, vcmrot.eigval) + oneT)), M)
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
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  fisher_B!(H, TwoVarCompModelRotate(vcm), vcobsrot)
end

function fisher_B{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  H = zeros(T, nmeanparams(vcm), nmeanparams(vcm))
  fisher_B!(H, TwoVarCompModelRotate(vcm), vcobsrot)
end

function fisher_B!{T <: AbstractFloat}(
  H::AbstractMatrix{T},
  vcm::VarianceComponentModel{T, 2},
  vcobs::VarianceComponentVariate{T, 2}
  )

  fisher_B!(H, TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs))
end

function fisher_B{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobs::VarianceComponentVariate{T, 2}
  )

  H = zeros(T, nmeanparams(vcm), nmeanparams(vcm))
  fisher_B!(H, TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs))
end

function fisher_B!{T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate},
  T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate}}(

  H::AbstractMatrix,
  vcm::T1,
  vcobs::Array{T2}
  )

  copy!(H, fisher_B(vcm, vcobs))
end

function fisher_B{T1 <: Union{VarianceComponentModel, TwoVarCompModelRotate},
  T2 <: Union{VarianceComponentVariate, TwoVarCompVariateRotate}}(

  vcm::T1,
  vcobs::Array{T2}
  )

  mapreduce(x -> fisher_B(vcm, x), +, vcobs)
end

#---------------------------------------------------------------------------#
# Update mean parameters by quadratic programming
#---------------------------------------------------------------------------#

"""
Compute the quadratic and linear parts of generalized least squares criterion.
"""
function suffstats{T <: AbstractFloat}(
  vcmrot::TwoVarCompModelRotate{T},
  vcdatarot::TwoVarCompVariateRotate{T}
  )

  oneT = one(T)
  wt = oneT ./ √(kron(vcmrot.eigval, vcdatarot.eigval) + oneT)
  Xy = [kron(vcmrot.eigvec', vcdatarot.Xrot) vec(vcdatarot.Yrot * vcmrot.eigvec)]
  scale!(wt, Xy)
  return At_mul_B(Xy, Xy)
end

function suffstats{T1 <: TwoVarCompModelRotate, T2 <: TwoVarCompVariateRotate}(
  vcmrot::T1,
  vcdatarot::Array{T2}
  )

  mapreduce(x -> suffstats(vcmrot, x), +, vcdatarot)
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
function update_meanparam!{T1 <: VarianceComponentModel,
  T2<:TwoVarCompVariateRotate}(
  vcm::T1,
  vcdatarot::Union{T2, Array{T2}},
  qpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver = MathProgBase.defaultQPsolver,
  Q::Matrix = Array{eltype(vcm)}(nmeanparams(vcm), nmeanparams(vcm)),
  c::Vector = Array{eltype(vcm)}(nmeanparams(vcm))
  )

  # quick return if there is no mean parameters
  isempty(vcm.B) && return vcm.B
  # rotate the model
  vcmrot = TwoVarCompModelRotate(vcm)
  # accumlate the quadratic and linear parts of QP
  Q = suffstats(vcmrot, vcdatarot)
  # quadratic programming
  qpsol = quadprog(-Q[1:end-1, end], Q[1:end-1, 1:end-1],
    vcm.A, vcm.sense, vcm.b, vcm.lb, vcm.ub, qpsolver)
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
type TwoVarCompOptProb{T1<:VarianceComponentModel,
  T2<:TwoVarCompVariateRotate} <: MathProgBase.AbstractNLPEvaluator
  # variance component model and data
  vcmodel::T1
  vcdatarot::Union{T2, Array{T2}}
  # QP solver for updating B given Σ
  qpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver
  # intermediate variables
  L::NTuple{2, Matrix} # Cholesky factors
  ∇Σ::Vector           # graident wrt (Σ1, Σ2)
  HΣ::Matrix           # Hessian wrt (Σ1, Σ2)
  HL::Matrix           # Hessian wrt (L1, L2)
  Q::Matrix            # pd-by-pd working matrix
  c::Vector            # pd working vector
end

"""
    TwoVarCompOptProb(vcm, vcobsrot)
    TwoVarCompOptProb(vcm, vcobsrot, qpsolver)

Constructor of [`TwoVarCompOptProb`](@ref) from [`VarianceComponentModel`](@ref) `vcm`,
[`TwoVarCompVariateRotate`](@ref) `vcobsrot`, and input quadratic programming
sovler `qpsolver`.
"""
function TwoVarCompOptProb{T1<:VarianceComponentModel, T2<:TwoVarCompVariateRotate}(
  vcm::T1,
  vcobsrot::Union{T2, Array{T2}},
  qpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver = MathProgBase.defaultQPsolver
  )

  T = eltype(vcm)
  d, pd = length(vcm), nmeanparams(vcm)
  # number of optimization parameters in variance
  nvar = nvarparams(vcm)
  # allocate intermediate variables
  L  = (zeros(T, d, d), zeros(T, d, d))
  ∇Σ = zeros(T, 2d^2) # graident wrt (Σ1, Σ2)
  HΣ = zeros(T, 2d^2, 2d^2) # Hessian wrt (Σ1, Σ2)
  HL = zeros(T, nvar, nvar) # Hessian wrt Ls
  Q  = zeros(T, pd, pd)
  c  = zeros(T, pd)
  TwoVarCompOptProb(vcm, vcobsrot, qpsolver, L, ∇Σ, HΣ, HL, Q, c)
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
  update_meanparam!(dd.vcmodel, dd.vcdatarot, dd.qpsolver, dd.Q, dd.c)
  # evaluate profile log-pdf
  logpdf(dd.vcmodel, dd.vcdatarot)
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
  update_meanparam!(dd.vcmodel, dd.vcdatarot, dd.qpsolver, dd.Q, dd.c)
  # gradient wrt (Σ[1], Σ[2])
  gradient!(dd.∇Σ, dd.vcmodel, dd.vcdatarot)
  # chain rule for gradient wrt Cholesky factor
  chol_gradient!(sub(grad_f, 1:nvarhalf), dd.∇Σ[1:d^2], dd.L[1])
  chol_gradient!(sub(grad_f, (nvarhalf+1):nvar),
    dd.∇Σ[(d^2+1):end], dd.L[2])
  # chain rule for exponential of diagonal entries
  for j in 1:d
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
  d = size(dd.L[1], 1)
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
  update_meanparam!(dd.vcmodel, dd.vcdatarot, dd.qpsolver, dd.Q, dd.c)
  # Hessian wrt (Σ1, Σ2)
  fisher!(dd.HΣ, dd.vcmodel, dd.vcdatarot)
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
  for j in 1:d
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
  update_meanparam!(vcmodel, vcdatarot, dd.qpsolver, dd.Q, dd.c)
  # update final objective value
  maxlogl = logpdf(vcmodel, vcdatarot)

  # standard errors
  Bcov = zeros(T, nmean, nmean)
  Bcov = inv(fisher_B(vcmodel, vcdatarot))
  Bse = similar(vcmodel.B)
  copy!(Bse, sqrt(diag(Bcov)))
  Σcov = inv(fisher(vcmodel, vcdatarot))
  Σse = (zeros(T, d, d), zeros(T, d, d))
  copy!(Σse[1], sqrt(diag(sub(Σcov, 1:d^2, 1:d^2))))
  copy!(Σse[2], sqrt(diag(sub(Σcov, d^2+1:2d^2, d^2+1:2d^2))))

  # output
  maxlogl, vcmodel, Σse, Σcov, Bse, Bcov
end # function mle_fs

#---------------------------------------------------------------------------#
# MM algorithm
#---------------------------------------------------------------------------#

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

  # initialize algorithm
  T = eltype(vcm)
  # n = no. observations, d = no. categories
  d, pd = length(vcm), nmeanparams(vcm)
  zeroT, oneT, halfT = zero(T), one(T), convert(T, 0.5)
  nmean = nmeanparams(vcm)
  # set up QP solver
  if qpsolver == :Ipopt
    qs = IpoptSolver(print_level = 0)
  elseif qpsolver == :Gurobi
    qs = GurobiSolver(OutputFlag = 0)
  elseif qpsolver == :Mosek
    qs = MosekSolver(MSK_IPAR_LOG = 0)
  end
  # initial log-likelihood
  vcmrot = TwoVarCompModelRotate(vcm)
  logl::T = logpdf(vcmrot, vcdatarot)
  if verbose
    println()
    println("     MM Algorithm")
    println("  Iter      Objective  ")
    println("--------  -------------")
    @printf("%8.d  %13.e\n", 0, logl)
  end
  # allocate intermediate variables
  if nmean > 0
    Q = zeros(T, pd, pd)
    c = zeros(T, pd)
  end
  Wt = oneT ./ √(vcdatarot.eigval * vcmrot.eigval' + oneT)
  res = residual(vcmrot, vcdatarot) .* Wt
  Whalf = zeros(T, n, d)
  dg = zeros(T, d)
  λj = zeroT

  # MM loop
  for iter = 1:maxiter

    # update Σ1
    @inbounds for j in 1:d
      λj = vcmrot.eigval[j]
      dg[j] = mapreduce(x -> x / (λj * x + oneT), +, vcdatarot.eigval)
      dg[j] = sqrt(dg[j])
    end
    Whalf = res .* Wt
    scale!(sqrt(vcdatarot.eigval), Whalf)
    scale!(Whalf, vcmrot.eigval .* dg)
    #W = sqrtm(Whalf' * Whalf) # produces imaginery eigenvalues due to precision
    #dg = 1.0 ./ dg
    #Σ[1] = scale(dg, inv(Φ))' * W * scale(dg, inv(Φ))
    # this approach is more numerical stable
    Whalfsvd = svdfact(Whalf)
    copy!(vcm.Σ[1], scale(sqrt(Whalfsvd[:S]),
      Whalfsvd[:Vt]) * scale(oneT ./ dg, inv(vcmrot.eigvec)))
    copy!(vcm.Σ[1], At_mul_B(vcm.Σ[1], vcm.Σ[1]))

    # update Σ2
    @inbounds for j = 1:d
      λj = vcmrot.eigval[j]
      #dg[j] = sqrt(sum(oneT ./ (vcmrot.eigval[j] * vcdatarot.eigval + oneT)))
      dg[j] = mapreduce(x -> oneT / (λj * x + oneT), +, vcdatarot.eigval)
      dg[j] = sqrt(dg[j])
    end
    Whalf = res .* Wt
    scale!(Whalf, dg)
    # W = sqrtm(Whalf' * Whalf)
    # dg = 1.0 ./ dg
    # Σ[2] = scale(dg, inv(Φ))' * W * scale(dg, inv(Φ))
    # this approach is more numerical stable
    Whalfsvd = svdfact(Whalf)
    copy!(vcm.Σ[2], scale(sqrt(Whalfsvd[:S]),
      Whalfsvd[:Vt]) * scale(oneT ./ dg, inv(vcmrot.eigvec)))
    copy!(vcm.Σ[2], At_mul_B(vcm.Σ[2], vcm.Σ[2]))
    # make sure the last varianec component is pos. def.
    ϵ = convert(T, 1e-8)
    clamp_diagonal!(vcm.Σ[2], ϵ, T(Inf))

    # update mean parameters
    if !isempty(vcm.B)
      update_meanparam!(vcm, vcdatarot, qs, Q, c)
    end
    vcmrot = TwoVarCompModelRotate(vcm)
    Wt = oneT ./ √(vcdatarot.eigval * vcmrot.eigval' + oneT)
    res = residual(vcmrot, vcdatarot) .* Wt

    # check convergence
    loglold = logl
    logl    = logpdf(vcmrot, vcdatarot)
    if verbose
      if (iter <= 10) || (iter > 10 && iter % 10 == 0)
        @printf("%8.d  %13.e\n", iter, logl)
      end
    end
    if abs(logl - loglold) < funtol * (abs(logl) + oneT)
      break
    end
  end
  if verbose; println(); end

  # standard errors
  Bcov = zeros(T, nmean, nmean)
  fisher_B!(Bcov, vcm, vcdatarot)
  Bcov = inv(Bcov)
  Bse = similar(vcm.B)
  copy!(Bse, sqrt(diag(Bcov)))
  Σcov = inv(fisher(vcm, vcdatarot))
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
function fit_mle!{T <: AbstractFloat}(
  vcmodel::VarianceComponentModel{T, 2},
  vcdata::VarianceComponentVariate{T, 2};
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
function fit_reml!{T <: AbstractFloat}(
  vcmodel::VarianceComponentModel{T, 2},
  vcdata::VarianceComponentVariate{T, 2};
  algo::Symbol = :FS
  )

  n, d = size(vcdata.Y, 1), size(vcdata.Y, 2)
  vcdatarot = TwoVarCompVariateRotate(vcdata)
  resrot = vcdatarot.Yrot - vcdatarot.Xrot * (vcdatarot.Xrot \ vcdatarot.Yrot)
  # use residuals as responses
  resdatarot = TwoVarCompVariateRotate(resrot, zeros(T, n, 0),
    vcdatarot.eigval, vcdatarot.logdetV2)
  resmodel = VarianceComponentModel(zeros(T, 0, d), vcmodel.Σ)
  if algo == :FS
    _, _, Σse, Σcov, = mle_fs!(resmodel, resdatarot)
  elseif algo == :MM
    _, _, Σse, Σcov, = mle_mm!(resmodel, resdatarot)
  end

  # estimate mean parameters from covariance estimate
  oneT = one(T)
  resmodelrot = TwoVarCompModelRotate(resmodel)
  wt = oneT ./ sqrt(kron(resdatarot.eigval, resmodelrot.eigval) + oneT)
  Xnew = kron(resmodelrot.eigvec', vcdatarot.Xrot)
  Ynew = vec(vcdatarot.Yrot * resmodelrot.eigvec)
  scale!(wt, Xnew)
  Ynew = Ynew .* wt
  copy!(vcmodel.B, Xnew \ Ynew)

  # standard errors and covariance of mean parameters
  nmean = nmeanparams(vcmodel)
  covmatrix = zeros(T, nmean + 2d^2, nmean + 2d^2)
  fisher!(covmatrix, vcmodel, vcdatarot)
  Bcov = inv(covmatrix[1:nmean, 1:nmean])
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

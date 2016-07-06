import Base.gradient

export heritability,
  logpdf, gradient!, gradient, fisher!, fisher,
  mle_fs!, mle_mm!,
  fit_mle!, fit_reml!,
  update_meanparam!

#---------------------------------------------------------------------------#
# Evaluate log-pdf
#---------------------------------------------------------------------------#

"""
Calculate log-pdf of a *rotated* two variance component instance under a
*rotated* two variance component model.
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
Calculate log-pdf of a *rotated* two variance component instance under an
*unrotated* two variance component model.
"""
function logpdf{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  logpdf(TwoVarCompModelRotate(vcm), vcobsrot)
end

"""
Calculate log-pdf of an *unrotated* two variance component instance under an
*unrotated* two variance component model.
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
# function logpdf(
#   vcm::Union{TwoVarCompModelRotate, VarianceComponentModel},
#   vcobs::Union{Array{TwoVarCompVariateRotate}, Array{VarianceComponentVariate}}
#   )
#
#   mapreduce(x -> logpdf(vcm, x), +, vcobs)
# end

function logpdf(
  vcm::TwoVarCompModelRotate,
  vcobs::Array{TwoVarCompVariateRotate}
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
- `vcmrot::TwoVarCompModelRotate`: *rotated* two variance component model.
- `vcobsrot::TwoVarCompVariateRotate`: *rotated* two variance component data.

# Output
- `∇`: gradient vector at `Σ = (Σ[1], Σ[2])`.

# TODO: optimize computation
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

function gradient{T <: AbstractFloat}(
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  d = length(vcmrot)
  ∇ = zeros(T, 2d^2)
  gradient!(∇, vcmrot, vcobsrot)
end

function gradient!{T <: AbstractFloat}(
  ∇::AbstractVector{T},
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::Array{TwoVarCompVariateRotate{T}}
  )

  fill!(∇, zero(T))
  tmp = copy(∇)
  @inbounds for i in eachindex(vcobsrot)
    gradient!(tmp, vcmrot, vcobsrot[i])
    ∇ += tmp
  end
  ∇
end

function gradient!{T <: AbstractFloat}(
  ∇::AbstractVector{T},
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  gradient!(∇, TwoVarCompModelRotate(vcm), vcobsrot)
end

function gradient!{T <: AbstractFloat}(
  ∇::AbstractVector{T},
  vcm::VarianceComponentModel{T, 2},
  vcobs::VarianceComponentVariate{T, 2}
  )

  gradient!(∇, TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs))
end

function gradient{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  gradient(TwoVarCompModelRotate(vcm), vcobsrot)
end

function gradient{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobs::VarianceComponentVariate{T, 2}
  )

  gradient(TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobs))
end

function gradient{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::Array{TwoVarCompVariateRotate{T}}
  )

  d = length(vcmrot)
  ∇ = zeros(T, 2d^2)
  gradient!(∇, vcmrot, vcobsrot)
end

#---------------------------------------------------------------------------#
# Evaluate Fisher information matrix
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
  C[:] = [mapreduce(
    x -> x^2 / (vcmrot.eigval[i] * x + oneT) / (vcmrot.eigval[j] * x + oneT),
    +, vcobsrot.eigval) for i=1:d, j=1:d]
  A_mul_Bt!(sub(H, 1:d^2, 1:d^2), scale(Φ2, vec(C)), Φ2)
  # (2, 1) block
  C[:] = [mapreduce(
    x -> x / (vcmrot.eigval[i] * x + oneT) / (vcmrot.eigval[j] * x + oneT), +,
    vcobsrot.eigval) for i=1:d, j=1:d]
  A_mul_Bt!(sub(H, (d^2+1):(2d^2), 1:d^2), scale(Φ2, vec(C)), Φ2)
  # d-by-d (2, 2) block
  C[:] = [mapreduce(
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

function fisher!{T, BT, YT, XT}(
  H::AbstractMatrix{T},
  vcmrot::TwoVarCompModelRotate{T, BT},
  vcobsrot::Array{TwoVarCompVariateRotate{T, YT, XT}}
  )

  fill!(H, zero(T))
  tmp = copy(H)
  for i in eachindex(vcobsrot)
    fisher!(tmp, vcmrot, vcobsrot)
    H += tmp
  end
  H
end

function fisher{T <: AbstractFloat}(
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::Array{TwoVarCompVariateRotate{T}}
  )

  d = length(vcmrot.eigval)
  H = zeros(T, 2d^2, 2d^2)
  fisher!(H, vcmrot, vcobsrot)
end

"""

    fisher_B!(H, vcmrot, vcobsrot)

Calculate Fisher information matrix of `B` and overwrite `H`
based on two variance component model `vcm` and *rotated* two
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

  H = zeros(T, nmeanparams(vcm), nmeanparams(vcmrot))
  fisher_B!(H, TwoVarCompModelRotate(vcm), vcobsrot)
end

#---------------------------------------------------------------------------#
# Update mean parameters by quadratic programming
#---------------------------------------------------------------------------#

"""
Update mean parameters `vcm.B` by quadratic programming.

## Input
- `vcm::VarianceComponentModel{T, 2}`: two variance component model
- `vcdatarot::TwoVarCompVariateRotate`: *rotated* two variance component data
- `qpsolver`: QP solver
- `Xwork`: working design matrix for updating `vcm.B`
- `ywork`: working responses for updating `vcm.B`
- `Q`: quadratic part in the QP
- `c`: linear part in the QP

## Output
- `vcm.B`: updated mean parameters
"""
function update_meanparam!{T <: AbstractFloat}(
    vcm::VarianceComponentModel{T, 2},
    vcdatarot::TwoVarCompVariateRotate{T},
    qpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver = MathProgBase.defaultQPsolver,
    Xwork::Matrix{T} = Array{T}(length(vcdatarot.Yrot), nmeanparams(vcm)),
    ywork::Vector{T} = Array{T}(length(vcdatarot.Yrot)),
    Q::Matrix{T} = Array{T}(nmeanparams(vcm), nmeanparams(vcm)),
    c::Vector{T} = Array{T}(nmeanparams(vcm))
    )

    # quick return if there is no mean parameters
    isempty(vcm.B) && return vcm.B
    # fill Xwork, ywork
    zeroT, oneT, infT = zero(T), one(T), convert(T, Inf)
    vcmrot = TwoVarCompModelRotate(vcm)
    wt = oneT ./ √(kron(vcmrot.eigval, vcdatarot.eigval) + oneT)
    fill!(Xwork, zeroT)
    kronaxpy!(vcmrot.eigvec', vcdatarot.Xrot, Xwork)
    copy!(ywork, vcdatarot.Yrot * vcmrot.eigvec)
    scale!(wt, Xwork)
    ywork = ywork .* wt
    # no constraints: quick return
    if size(vcm.A, 1) == 0 && all(vcm.lb .== -infT) && all(vcm.ub .== infT)
      copy!(vcm.B, Xwork \ ywork)
      return vcm.B
    end
    # constraints: quadratic programming
    At_mul_B!(Q, Xwork, Xwork)
    At_mul_B!(c, Xwork, ywork)
    scale!(c, -oneT)
    qpsol = quadprog(c, Q, vcm.A, vcm.sense, vcm.b, vcm.lb, vcm.ub, qpsolver)
    if qpsol.status ≠ :Optimal
      println("Error in quadratic programming $(qpsol.status)")
    end
    copy!(vcm.B, qpsol.sol)
    return vcm.B
end

#---------------------------------------------------------------------------#
# Fisher scoring algorithm
#---------------------------------------------------------------------------#

type TwoVarCompOptProb{T <: AbstractFloat} <: MathProgBase.AbstractNLPEvaluator
  vcmodel::VarianceComponentModel{T, 2}
  vcdatarot::Union{TwoVarCompVariateRotate{T}, Array{TwoVarCompVariateRotate{T}}}
  # QP solver for updating B given Σ
  qpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver
  # intermediate variables
  L::NTuple{2, Matrix{T}} # Cholesky factors
  ∇Σ::Vector{T} # graident wrt (Σ1, Σ2)
  HΣ::Matrix{T} # Hessian wrt (Σ1, Σ2)
  HL::Matrix{T} # Hessian wrt (L1, L2)
  Xwork::Matrix{T} # nd-by-pd working matrix
  ywork::Vector{T} # nd working vector
  Q::Matrix{T}     # pd-by-pd working matrix
  c::Vector{T}     # pd working vector
end

function TwoVarCompOptProb{T}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::Union{TwoVarCompVariateRotate{T}, Array{TwoVarCompVariateRotate{T}}},
  qpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver = MathProgBase.defaultQPsolver
  )

  n, d, p = size(vcobsrot.Yrot, 1), size(vcobsrot.Yrot, 2), size(vcobsrot.Xrot, 2)
  # number of optimization parameters in variance
  nvar = nvarparams(vcm)
  # allocate intermediate variables
  zeroT = convert(T, 0)
  L = (zeros(T, d, d), zeros(T, d, d))
  ∇Σ = zeros(T, 2d^2) # graident wrt (Σ1, Σ2)
  HΣ = zeros(T, 2d^2, 2d^2) # Hessian wrt (Σ1, Σ2)
  HL = zeros(T, nvar, nvar) # Hessian wrt Ls
  Xwork = zeros(T, n * d, p * d)
  ywork = zeros(T, n * d)
  Q = zeros(T, p * d, p * d)
  c = zeros(T, p * d)
  TwoVarCompOptProb{T}(vcm, vcobsrot, qpsolver, L, ∇Σ, HΣ, HL, Xwork, ywork, Q, c)
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
  d = size(dd.L[1], 1)
  nvar = nvarparams(dd.vcmodel)
  nvarhalf = div(nvar, 2)
  # variance parameter
  dd.L[1][trilind(dd.L[1])] = x[1:nvarhalf]
  dd.L[2][trilind(dd.L[2])] = x[(nvarhalf+1):end]
  A_mul_Bt!(dd.vcmodel.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(dd.vcmodel.Σ[2], dd.L[2], dd.L[2])
  # update mean parameters
  update_meanparam!(dd.vcmodel, dd.vcdatarot,
    dd.qpsolver, dd.Xwork, dd.ywork, dd.Q, dd.c)
  # evaluate profile log-pdf
  logpdf(dd.vcmodel, dd.vcdatarot)
end # function MathProgBase.eval_f

function MathProgBase.eval_grad_f{T}(
  dd::TwoVarCompOptProb,
  grad_f::Vector{T},
  x::Vector{T}
  )

  d = size(dd.L[1], 1)
  nvar = nvarparams(dd.vcmodel)
  nvarhalf = div(nvar, 2)
  # variance parameter
  dd.L[1][trilind(dd.L[1])] = x[1:nvarhalf]
  dd.L[2][trilind(dd.L[2])] = x[(nvarhalf+1):end]
  A_mul_Bt!(dd.vcmodel.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(dd.vcmodel.Σ[2], dd.L[2], dd.L[2])
  # update mean parameters
  update_meanparam!(dd.vcmodel, dd.vcdatarot,
    dd.qpsolver, dd.Xwork, dd.ywork, dd.Q, dd.c)
  # gradient wrt (Σ[1], Σ[2])
  gradient!(dd.∇Σ, dd.vcmodel, dd.vcdatarot)
  # chain rule for gradient wrt Cholesky factor
  chol_gradient!(sub(grad_f, 1:nvarhalf), dd.∇Σ[1:d^2], dd.L[1])
  chol_gradient!(sub(grad_f, (nvarhalf+1):nvar),
    dd.∇Σ[(d^2+1):end], dd.L[2])
end # function MathProgBase.eval_grad_f

function MathProgBase.hesslag_structure(dd::TwoVarCompOptProb)
  d = size(dd.L[1], 1)
  nvar = nvarparams(dd.vcmodel)
  # linear indices for variance parameters
  ind2sub((nvar, nvar), trilind(nvar))
end # function MathProgBase.hesslag_structure

function MathProgBase.eval_hesslag{T}(dd::TwoVarCompOptProb, H::Vector{T},
  x::Vector{T}, σ::T, μ::Vector{T})

  d = size(dd.L[1], 1)
  nvar = nvarparams(dd.vcmodel)
  nvarhalf = div(nvar, 2)
  # variance parameter
  dd.L[1][trilind(dd.L[1])] = x[1:nvarhalf]
  dd.L[2][trilind(dd.L[2])] = x[(nvarhalf+1):end]
  A_mul_Bt!(dd.vcmodel.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(dd.vcmodel.Σ[2], dd.L[2], dd.L[2])
  # update mean parameters
  update_meanparam!(dd.vcmodel, dd.vcdatarot,
    dd.qpsolver, dd.Xwork, dd.ywork, dd.Q, dd.c)
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
  # output
  copy!(H, vech(dd.HL))
  scale!(H, -σ)
end

"""
    mle_fs!(vcmodel, vcdatarot; maxiter, solver, qpsolver, verbose)

Find MLE by Fisher scoring algorithm.
"""
function mle_fs!{T}(
  vcmodel::VarianceComponentModel{T, 2},
  vcdatarot::TwoVarCompVariateRotate{T};
  maxiter::Integer = 1000,
  solver::Symbol = :Ipopt,
  qpsolver::Symbol = :Ipopt,
  verbose::Bool = true,
  )

  n, d = size(vcdatarot.Yrot, 1), size(vcdatarot.Yrot, 2)
  nd = n * d
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
      KTR_PARAM_ALG = 1,
      KTR_PARAM_OUTLEV = verbose? 2 : 0,
      #KTR_PARAM_GRADOPT = 1,
      #KTR_PARAM_HESSOPT = 1,
      #KTR_PARAM_DERIVCHECK = 2
      #KTR_PARAM_TUNER = 1,
      #KTR_PARAM_MAXTIMECPU = 5.0,
      )
  end
  m = MathProgBase.NonlinearModel(solver)
  # lower and upper bounds for variance parameter
  lb = zeros(T, nvar)
  fill!(lb, convert(T, -Inf))
  for j in 1:d
    # linear index of diagonal entries of L1
    idx = 1 + (j - 1) * d - div((j - 1) * (j - 2), 2)
    lb[idx] = zeroT
    # linear index of diagonal entries of L2
    idx += binomial(d + 1, 2)
    lb[idx] = convert(T, 1e-4) # make sure last variance component is pos. def.
  end
  ub = similar(lb)
  fill!(ub, convert(T, Inf))
  MathProgBase.loadproblem!(m, nvar, 0, lb, ub, T[], T[], :Max, dd)
  # start point
  x0 = [vech(cholfact(vcmodel.Σ[1], :L, Val{true})[:L].data);
        vech(cholfact(vcmodel.Σ[2], :L, Val{true})[:L].data)]
  MathProgBase.setwarmstart!(m, x0)
  # optimize
  MathProgBase.optimize!(m)
  stat = MathProgBase.status(m)
  x = MathProgBase.getsolution(m)
  # retrieve variance component parameters
  dd.L[1][Ltrilind] = x[1:nvarhalf]
  dd.L[2][Ltrilind] = x[nvarhalf+1:end]
  A_mul_Bt!(vcmodel.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(vcmodel.Σ[2], dd.L[2], dd.L[2])
  # update mean parameters
  update_meanparam!(vcmodel, vcdatarot, dd.qpsolver,
    dd.Xwork, dd.ywork, dd.Q, dd.c)
  # update final objective value
  maxlogl = logpdf(vcmodel, vcdatarot)

  # standard errors
  Bcov = zeros(T, nmean, nmean)
  fisher_B!(Bcov, vcmodel, vcdatarot)
  Bcov = inv(Bcov)
  Bse = similar(vcmodel.B)
  copy!(Bse, sqrt(diag(Bcov)))
  Σcov = zeros(T, 2d^2, 2d^2)
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
  mle_mm!(Yrot, ev, loglconst; Σ0, maxiter, verbose)

Fit variance component model using minorization-maximization algorithm. Data
`vec(Y)` is assumed to be normal with mean zero and covariance
`Σ[1]⊗V[1] + Σ[2]⊗V[2]`.

# Input
- `Yrot`: rotated responses `U'*Y`, where `(ev,U) = eig(V[1],V[2])`.
- `ev`: eigenvalues from `(ev,U) = eig(V[1],V[2])`.
- `loglconst`: constant `n*d*log(2π)+d*logdet(V2)` in 2.0log-likelihood.

# Keyword arguments
- `Σ0=(Σ0[1], Σ0[2])`: starting value for variance component parameters.
- `maxiter`: maximum number of iterations for nonlinear programming solver.
- `verbose`: logging information.

# Output
- `logl`: log-likelihood at `Σ=(Σ[1],Σ[2])`.
- `Σ=(Σ[1],Σ[2])`: variance component estimates.
- `Σse=(Σse[1],Σse[2])`: standard errors of variance component estimates.
- `Σcov`: `2d^2 x 2d^2` covariance matrix of variance component estimates.

# Reference
- H. Zhou, L. Hu, J. Zhou, and K. Lange (2015)
  MM algorithms for variance components models.
  [http://arxiv.org/abs/1509.07426](http://arxiv.org/abs/1509.07426)
"""
function mle_mm!{T, BT, ΣT, YT, XT}(
  vcm::VarianceComponentModel{T, 2, BT, ΣT},
  vcdatarot::TwoVarCompVariateRotate{T, YT, XT};
  maxiter::Integer = 10000,
  funtol::T = convert(T, 1e-8),
  qpsolver::Symbol = :Ipopt,
  verbose::Bool = true,
  )

  # initialize algorithm
  # n = no. observations, d = no. categories
  n, d, p = size(vcdatarot.Yrot, 1), size(vcdatarot.Yrot, 2), size(vcm.B, 1)
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
    Xnew = zeros(T, n * d, p * d)
    Ynew = zeros(T, n * d)
    Q = zeros(T, p * d, p * d)
    c = zeros(T, p * d)
  end
  Wt = oneT ./ sqrt(vcdatarot.eigval * vcmrot.eigval' + oneT)
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

    # update mean parameters
    if !isempty(vcm.B)
      update_meanparam!(vcm, vcdatarot, qs, Xnew, Ynew, Q, c)
    end
    vcmrot = TwoVarCompModelRotate(vcm)
    Wt = oneT ./ sqrt(vcdatarot.eigval * vcmrot.eigval' + oneT)
    res = residual(vcmrot, vcdatarot) .* Wt

    # check convergence
    loglold = logl
    logl = logpdf(vcmrot, vcdatarot)
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
  Σcov = zeros(T, 2d^2, 2d^2)
  Σse = (zeros(T, d, d), zeros(T, d, d))
  copy!(Σse[1], sqrt(diag(sub(Σcov, 1:d^2, 1:d^2))))
  copy!(Σse[2], sqrt(diag(sub(Σcov, d^2+1:2d^2, d^2+1:2d^2))))

  # output
  logl, vcm, Σse, Σcov, Bse, Bcov
end # function mle_mm

#---------------------------------------------------------------------------#
# Estimation gateway
#---------------------------------------------------------------------------#

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

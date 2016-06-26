import Base.gradient

export heritability,
  logpdf, gradient!, gradient, fisher!, fisher,
  mle_fs!, mle_mm!

#---------------------------------------------------------------------------#
# Evaluate log-pdf
#---------------------------------------------------------------------------#

"""
Calculate log-pdf of a *rotated* two variance component instance under a
*rotated* two variance component model.
"""
function logpdf(
  vcmrot::TwoVarCompModelRotate,
  vcobsrot::TwoVarCompVariateRotate)

  n, d = size(vcobsrot.Yrot, 1), size(vcobsrot.Yrot, 2)
  nd = n * d
  T = eltype(vcmrot)
  zeroT, oneT = zero(T), one(T)
  res = residual(vcmrot, vcobsrot)
  # evaluate 2(log-likehood)
  objval::T = - nd * log(2π) - d * vcobsrot.logdetV2 - n * vcmrot.logdetΣ2
  tmp = zeroT
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
function logpdf{T}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::TwoVarCompVariateRotate
  )

  logpdf(TwoVarCompModelRotate(vcm), vcobsrot)
end

"""
Calculate log-pdf of an *unrotated* two variance component instance under an
*unrotated* two variance component model.
"""
function logpdf{T}(
  vcm::VarianceComponentModel{T, 2},
  vcobsrot::VarianceComponentVariate{T, 2}
  )

  logpdf(TwoVarCompModelRotate(vcm), TwoVarCompVariateRotate(vcobsrot))
end

"""
Calculate log-pdf of an array of two variance component instances
under a two variance component model.
"""
function logpdf{T}(
  vcm::Union{TwoVarCompModelRotate{T}, VarianceComponentModel{T, 2}},
  vcobs::Union{Array{TwoVarCompVariateRotate{T}}, Array{VarianceComponentVariate{T, 2}}}
  )

  mapreduce(x -> logpdf(vcm, x), +, vcobs)
end

#---------------------------------------------------------------------------#
# Evaluate gradient
#---------------------------------------------------------------------------#

"""

    gradient!(∇, vcmrot, vcobsrot)

Evaluate gradient at `Σ = (Σ[1], Σ[2])` and overwrite `∇`.

# Input
- `∇`: gradient vector.
- `vcmrot`: *rotated* two variance component model.
- `vcobsrot`: *rotated* two variance component data instance.

# Output
- `∇`: gradient vector at `Σ = (Σ[1], Σ[2])`.
"""
function gradient!{T <: AbstractFloat}(
  ∇::AbstractVector{T},
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  n, d = size(vcobsrot.Yrot, 1), size(vcobsrot.Yrot, 2)
  zeroT, oneT = zero(T), one(T)
  res = residual(vcmrot, vcobsrot)
  # evaluate gradient
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
  ∇[d^2+1:2d^2] = N2[:]
  scale!(∇, convert(T, 0.5))
end # function gradient!

function gradient{T <: AbstractFloat}(
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  d = length(vcmrot.eigval)
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

  d = length(vcmrot.eigval)
  ∇ = zeros(T, 2d^2)
  gradient!(∇, vcmrot, vcobsrot)
end

#---------------------------------------------------------------------------#
# Evaluate Fisher information matrix
#---------------------------------------------------------------------------#

"""

    fisher!(H, Σ, ev)

Calculate Fisher information matrix at `Σ = (Σ[1], Σ[2])` and overwrite `H`,
under the model `vec(Y)` is normal with mean zero and covariance
`Σ[1]⊗V[1] + Σ[2]⊗V[2]`.

# Input
- `H`: Hessian matrix.
- `Σ = (Σ[1], Σ[2])`: variance component parameters.
- `ev`: eigenvalues from `(λ, U) = eig(V1, V2)`.

# Output
- `H`: Fisher information matrix at `Σ = (Σ[1], Σ[2])`.
"""
function fisher!{T <: AbstractFloat}(
  H::AbstractMatrix{T},
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::TwoVarCompVariateRotate{T}
  )

  d = length(vcmrot.eigval)
  zeroT, oneT = zero(T), one(T)
  # evaluate Hessian
  C = zeros(T, d, d)
  Φ2 = kron(vcmrot.eigvec, vcmrot.eigvec)
  # (1, 1) block
  for l2 in 1:d, l1 in l2:d# only the lower triangular part
    C[l1, l2] = sum(vcobsrot.eigval .* vcobsrot.eigval ./
      (vcmrot.eigval[l1] * vcobsrot.eigval + oneT) ./
      (vcmrot.eigval[l2] * vcobsrot.eigval + oneT))
  end
  LinAlg.copytri!(C, 'L') # copy to upper triangular part
  A_mul_Bt!(sub(H, 1:d^2, 1:d^2), scale(Φ2, vec(C)), Φ2)
  # (2, 1) block
  for l2 in 1:d, l1 in l2:d
    C[l1, l2] = sum(vcobsrot.eigval ./
      (vcmrot.eigval[l1] * vcobsrot.eigval + oneT) ./
      (vcmrot.eigval[l2] * vcobsrot.eigval + oneT))
  end
  LinAlg.copytri!(C, 'L') # copy to upper triangular part
  A_mul_Bt!(sub(H, d^2+1:2d^2, 1:d^2), scale(Φ2, vec(C)), Φ2)
  # d-by-d (2, 2) block
  for l2 in 1:d, l1 in l2:d
    C[l1, l2] = sum(oneT ./ (vcmrot.eigval[l1] * vcobsrot.eigval + oneT) ./
      (vcmrot.eigval[l2] * vcobsrot.eigval + oneT))
  end
  LinAlg.copytri!(C, 'L') # copy to upper triangular part
  A_mul_Bt!(sub(H, d^2+1:2d^2, d^2+1:2d^2), scale(Φ2, vec(C)), Φ2)
  # make sure it's Hessian of the *negative* log-likehood
  LinAlg.copytri!(H, 'L') # copy to upper triangular part
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

function fisher!{T <: AbstractFloat}(
  H::AbstractMatrix{T},
  vcmrot::TwoVarCompModelRotate{T},
  vcobsrot::Array{TwoVarCompVariateRotate{T}}
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

#---------------------------------------------------------------------------#
# Fisher scoring algorithm
#---------------------------------------------------------------------------#

type TwoVarCompOptProb{T <: AbstractFloat} <: MathProgBase.AbstractNLPEvaluator
  vcmodel::VarianceComponentModel{T, 2}
  vcmodelrot::TwoVarCompModelRotate{T}
  vcdatarot::Union{TwoVarCompVariateRotate{T}, Array{TwoVarCompVariateRotate{T}}}
  L::NTuple{2, Matrix{T}}
  ∇Σ::Vector{T} # graident wrt Σs
  HΣ::Matrix{T} # Hessian wrt Σs
  HL::Matrix{T} # Hessian wrt Ls
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
  nparam = d * (d + 1)
  nparamhalf = div(nparam, 2)
  dd.L[1][trilind(dd.L[1])] = x[1:nparamhalf]
  dd.L[2][trilind(dd.L[2])] = x[nparamhalf+1:end]
  A_mul_Bt!(dd.vcmodel.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(dd.vcmodel.Σ[2], dd.L[2], dd.L[2])
  dd.vcmodelrot = TwoVarCompModelRotate(dd.vcmodel)
  logpdf(dd.vcmodelrot, dd.vcdatarot)
end # function MathProgBase.eval_f

function MathProgBase.eval_grad_f{T}(
  dd::TwoVarCompOptProb,
  grad_f::Vector{T},
  x::Vector{T}
  )

  d = size(dd.L[1], 1)
  nparam = d * (d + 1)
  nparamhalf = div(nparam, 2)
  dd.L[1][trilind(dd.L[1])] = x[1:nparamhalf]
  dd.L[2][trilind(dd.L[2])] = x[nparamhalf+1:end]
  A_mul_Bt!(dd.vcmodel.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(dd.vcmodel.Σ[2], dd.L[2], dd.L[2])
  # gradient wrt Σ
  dd.vcmodelrot = TwoVarCompModelRotate(dd.vcmodel)
  gradient!(dd.∇Σ, dd.vcmodel, dd.vcdatarot)
  # chain rule for gradient wrt Cholesky factor
  chol_gradient!(sub(grad_f, 1:nparamhalf),
    dd.∇Σ[1:d^2], dd.L[1])
  chol_gradient!(sub(grad_f, nparamhalf+1:nparam),
    dd.∇Σ[d^2+1:end], dd.L[2])
end # function MathProgBase.eval_grad_f

function MathProgBase.hesslag_structure(dd::TwoVarCompOptProb)
  d = size(dd.L[1], 1)
  nparam = d * (d + 1)
  ind2sub((nparam, nparam), trilind(nparam))
end # function MathProgBase.hesslag_structure

function MathProgBase.eval_hesslag{T}(dd::TwoVarCompOptProb, H::Vector{T},
  x::Vector{T}, σ::T, μ::Vector{T})
  d = size(dd.L[1], 1)
  nparam = d * (d + 1)
  nparamhalf = div(nparam, 2)
  dd.L[1][trilind(dd.L[1])] = x[1:nparamhalf]
  dd.L[2][trilind(dd.L[2])] = x[nparamhalf+1:end]
  A_mul_Bt!(dd.vcmodel.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(dd.vcmodel.Σ[2], dd.L[2], dd.L[2])
  dd.vcmodelrot = TwoVarCompModelRotate(dd.vcmodel)
  fisher!(dd.HΣ, dd.vcmodelrot, dd.vcdatarot)
  # chain rule for Hessian wrt Cholesky factor
  # only the lower left triangle
  # (1, 1) block
  chol_gradient!(sub(dd.HL, 1:nparamhalf, 1:nparamhalf),
    chol_gradient(dd.HΣ[1:d^2, 1:d^2], dd.L[1])', dd.L[1])
  # (2, 1) block
  chol_gradient!(sub(dd.HL, nparamhalf+1:nparam, 1:nparamhalf),
    chol_gradient(dd.HΣ[d^2+1:2d^2, 1:d^2], dd.L[1])', dd.L[2])
  # (2, 2) block
  chol_gradient!(sub(dd.HL, nparamhalf+1:nparam, nparamhalf+1:nparam),
    chol_gradient(dd.HΣ[d^2+1:2d^2, d^2+1:2d^2], dd.L[2])', dd.L[2])
  # output
  scale!(dd.HL, -σ)
  copy!(H, vech(dd.HL))
end

function mle_fs!{T <: AbstractFloat}(
  vcmodel::VarianceComponentModel{T, 2},
  vcdatarot::TwoVarCompVariateRotate{T};
  maxiter::Integer = 1000,
  solver::Symbol = :Ipopt,
  verbose::Bool = true
  )

  n, d = size(vcdatarot.Yrot, 1), size(vcdatarot.Yrot, 2)
  nd = n * d
  Ltrilind = trilind(d, d)
  # number of optimization parameters in Cholesky factors
  nparam = d * (d + 1)
  nparamhalf = div(nparam, 2)

  # pre-allocate variables for optimization
  zeroT = convert(T, 0)
  L = (zeros(T, d, d), zeros(T, d, d))
  ∇Σ = zeros(T, 2d^2) # graident wrt Σs
  HΣ = zeros(T, 2d^2, 2d^2) # Hessian wrt Σs
  HL = zeros(T, nparam, nparam) # Hessian wrt Ls
  # data for the optimization problem
  dd = TwoVarCompOptProb(vcmodel, TwoVarCompModelRotate(vcmodel),
    vcdatarot, L, ∇Σ, HΣ, HL)

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
      )
  end
  m = MathProgBase.NonlinearModel(solver)
  # lower and upper bounds
  lb = zeros(T, nparam)
  fill!(lb, convert(T, -Inf))
  for j in 1:d
    idx = 1 + (j - 1) * d - div((j - 1) * (j - 2), 2)
    lb[idx] = zeroT
    idx += div(d * (d + 1), 2)
    lb[idx] = convert(T, 1e-4) # make sure last variance component is pos. def.
  end
  ub = similar(lb)
  fill!(ub, convert(T, Inf))
  MathProgBase.loadproblem!(m, nparam, 0, lb, ub, T[], T[], :Max, dd)
  # start point
  x0 = [vech(chol(vcmodel.Σ[1], Val{:L}).data);
        vech(chol(vcmodel.Σ[2], Val{:L}).data)]
  MathProgBase.setwarmstart!(m, x0)
  # convergence criteria
  #xtol_rel!(opt, 1e-8)
  #ftol_rel!(opt, 1e-8)
  #maxtime!(opt, 60)
  # optimize
  MathProgBase.optimize!(m)
  stat = MathProgBase.status(m)
  x = MathProgBase.getsolution(m)
  maxlogl = MathProgBase.getobjval(m)
  # retrieve result
  L[1][Ltrilind] = x[1:nparamhalf]
  L[2][Ltrilind] = x[nparamhalf+1:end]
  A_mul_Bt!(vcmodel.Σ[1], L[1], L[1])
  A_mul_Bt!(vcmodel.Σ[2], L[2], L[2])

  # standard errors
  Σcov = zeros(T, 2d^2, 2d^2)
  fisher!(Σcov, vcmodel, vcdatarot)
  Σcov = inv(Σcov)
  Σse = deepcopy(vcmodel.Σ)
  copy!(Σse[1], sqrt(diag(sub(Σcov, 1:d^2, 1:d^2))))
  copy!(Σse[2], sqrt(diag(sub(Σcov, d^2+1:2d^2, d^2+1:2d^2))))

  # output
  maxlogl, Σse, Σcov
end # function mle_fs

#---------------------------------------------------------------------------#
# MM algorithm
#---------------------------------------------------------------------------#

"""
  reml_mm(Yrot, ev, loglconst; Σ0, maxiter, verbose)

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
function mle_mm!{T <: AbstractFloat}(
  vcm::VarianceComponentModel{T, 2},
  vcdatarot::Union{TwoVarCompVariateRotate{T}, Array{TwoVarCompVariateRotate{T}}};
  maxiter::Integer = 10000,
  funtol::T = convert(T, 1e-8),
  verbose::Bool = true)

  # initialize algorithm
  # n = no. observations, d = no. categories
  n, d = size(vcdatarot.Yrot, 1), size(vcdatarot.Yrot, 2)
  nd = n * d
  zeroT, oneT, halfT = zero(T), one(T), convert(T, 0.5)
  # update generalized eigen-decomposition
  vcmrot = TwoVarCompModelRotate(vcm)
  Wt = oneT ./ sqrt(vcdatarot.eigval * vcmrot.eigval' + oneT)
  res = (vcdatarot.Yrot * vcmrot.eigvec) .* Wt
  logl = sum(logpdf(vcmrot, vcdatarot))
  if verbose
    println()
    println("     MM Algorithm")
    println("  Iter      Objective  ")
    println("--------  -------------")
    @printf("%8.d  %13.e\n", 0, logl)
  end

  # MM loop
  Whalf = zeros(T, n, d)
  dg = zeros(T, d)
  for iter = 1:maxiter
    # update Σ1
    for j = 1:d
      @inbounds dg[j] = sqrt(sum(vcdatarot.eigval ./
        (vcmrot.eigval[j] * vcdatarot.eigval + oneT)))
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
    At_mul_B!(vcm.Σ[1], vcm.Σ[1], vcm.Σ[1])
    # update Σ2
    @inbounds for j = 1:d
      dg[j] = sqrt(sum(oneT ./ (vcmrot.eigval[j] * vcdatarot.eigval + oneT)))
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
    At_mul_B!(vcm.Σ[2], vcm.Σ[2], vcm.Σ[2])

    # update generalized eigen-decomposition
    vcmrot = TwoVarCompModelRotate(vcm)
    Wt = oneT ./ sqrt(vcdatarot.eigval * vcmrot.eigval' + oneT)
    res = (vcdatarot.Yrot * vcmrot.eigvec) .* Wt

    # check convergence
    loglold = logl
    logl = sum(logpdf(vcmrot, vcdatarot))
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
  Σcov = zeros(T, 2d^2, 2d^2)
  fisher!(Σcov, vcm, vcdatarot)
  Σcov = inv(Σcov)
  Σse = deepcopy(vcm.Σ)
  copy!(Σse[1], sqrt(diag(sub(Σcov, 1:d^2, 1:d^2))))
  copy!(Σse[2], sqrt(diag(sub(Σcov, d^2+1:2d^2, d^2+1:2d^2))))

  # output
  logl, Σse, Σcov
end # function mle_mm

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

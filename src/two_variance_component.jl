export reml_objval,
  reml_grad, reml_grad!,
  reml_fisher, reml_fisher!,
  reml_eig, reml_fs, reml_mm,
  heritability,
  logpdf

"""

    reml_objval(Σ, Yrot, ev, loglconst)

Evaluate log-likehood at `Σ = (Σ[1], Σ[2])` of the model `vec(Y)` is
normal with mean zero and covariance `Σ[1]⊗V[1] + Σ[2]⊗V[2]`.

# Input
- `Σ = (Σ[1], Σ[2])`: variance component parameters.
- `Yrot`: rotated responses `U' * Y`, where `(ev, U) = eig(V[1], V[2])`.
- `ev`: eigenvalues from `(λ, U) = eig(V[1], V[2])`.
- `loglconst`: constant `n*d*log(2π) + d*logdet(V2)` in 2.0log-likelihood.

# Output
- log-likelihood at `Σ = (Σ[1], Σ[2])`.
"""
function reml_objval{T <: AbstractFloat}(
  Σ::Union{Vector{Matrix{T}}, Tuple{Matrix{T}, Matrix{T}}},
  Yrot::AbstractVecOrMat{T},
  ev::Vector{T},
  loglconst::T
  )

  n, d = size(Yrot, 1), size(Yrot, 2)
  zeroT, oneT = zero(T), one(T)
  # generalized eigenvalue decomposition of (Γ1, Γ2) = (L1L1t, L2L2t)
  # normalize columns of Φ to be orthonormal with respect to Σ[2]
  λ, Φ = eig(Σ[1], Σ[2])
  map!(x -> max(x, zeroT), λ) # correct negative eigenvalues due to roundoff
  #scale!(Φ, 1.0 ./ sqrt(diag(Φ' * Σ[2] * Φ)))
  # rotate the columns of Yrot
  Ynew = Yrot * Φ
  # evaluate 2.0(log-likehood)
  objval = loglconst - n * logdet(Σ[2])
  tmp = zeroT
  @inbounds for j in 1:d
    λj = λ[j]
    @simd for i in 1:n
      tmp = oneT / (ev[i] * λj + oneT)
      objval += log(tmp) - tmp * Ynew[i, j]^2
    end
  end
  objval *= convert(T, 0.5)
end # function reml_objval

"""

    reml_grad!(∇, Σ, Yrot, ev)

Evaluate gradient at `Σ = (Σ[1], Σ[2])` and overwrite `∇`, under the model
`vec(Y)` is normal with mean zero and covariance `Σ[1]⊗V[1] + Σ[2]⊗V[2]`.

# Input
- `∇`: gradient vector.
- `Σ = (Σ[1], Σ[2])`: variance component parameters.
- `Yrot`: rotated responses `U' * Y`, where `(ev, U) = eig(V1, V2)`.
- `ev`: eigenvalues from `(λ, U) = eig(V1, V2)`.

# Output
- `∇`: gradient vector at `Σ = (Σ[1], Σ[2])`.
"""
function reml_grad!{T <: AbstractFloat}(
  ∇::Vector{T},
  Σ::Union{Vector{Matrix{T}}, Tuple{Matrix{T}, Matrix{T}}},
  Yrot::AbstractVecOrMat{T},
  ev::Vector{T}
  )

  n, d = size(Yrot, 1), size(Yrot, 2)
  zeroT, oneT = zero(T), one(T)
  # generalized eigenvalue decomposition of (Γ1, Γ2) = (L1L1t, L2L2t)
  # normalize columns of Φ to be orthonormal with respect to Σ[2]
  λ, Φ = eig(Σ[1], Σ[2])
  map!(x -> max(x, zeroT), λ) # correct negative eigenvalues due to roundoff
  #scale!(Φ, 1.0 ./ sqrt(diag(Φ' * Σ[2] * Φ)))
  # rotate the columns of Yrot
  Ynew = Yrot * Φ
  # evaluate gradient
  m1diag = zeros(T, d)
  m2diag = zeros(T, d)
  tmp = zeroT
  @inbounds for j in 1:d
    @simd for i in 1:n
      tmp = oneT / (ev[i] * λ[j] + oneT)
      Ynew[i, j] *= tmp
      m1diag[j] += ev[i] * tmp
      m2diag[j] += tmp
    end
  end
  N2 = At_mul_B(Ynew, Ynew)
  scale!(sqrt(ev), Ynew)
  N1 = At_mul_B(Ynew, Ynew)
  @inbounds for j in 1:d
    N1[j, j] -= m1diag[j]
    N2[j, j] -= m2diag[j]
  end
  N1 = Φ * N1 * Φ'
  N2 = Φ * N2 * Φ'
  ∇[1:d^2] = N1[:]
  ∇[d^2+1:2d^2] = N2[:]
  scale!(∇, convert(T, 0.5))
end # function reml_grad!

function reml_grad{T <: AbstractFloat}(
  Σ::Union{Vector{Matrix{T}}, Tuple{Matrix{T}, Matrix{T}}},
  Yrot::AbstractVecOrMat{T},
  ev::Vector{T}
  )

  d = size(Σ[1], 1)
  ∇ = zeros(T, 2d^2)
  reml_grad!(∇, Σ, Yrot, ev)
  ∇
end # function reml_grad

"""

    reml_fisher!(H, Σ, ev)

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
function reml_fisher!{T <: AbstractFloat}(
  H::Matrix{T},
  Σ::Union{Vector{Matrix{T}}, Tuple{Matrix{T}, Matrix{T}}},
  ev::Vector{T}
  )

  d = size(Σ[1], 1)
  zeroT, oneT = zero(T), one(T)
  # generalized eigenvalue decomposition of (Γ1, Γ2) = (L1L1t, L2L2t)
  # normalize columns of Φ to be orthonormal with respect to Σ[2]
  λ, Φ = eig(Σ[1], Σ[2])
  #scale!(Φ, 1.0 ./ sqrt(diag(Φ' * Σ[2] * Φ)))
  map!(x -> max(x, zeroT), λ) # correct negative eigenvalues due to roundoff
  # evaluate Hessian
  C = zeros(T, d, d)
  Φ2 = kron(Φ, Φ)
  # (1, 1) block
  for l2 in 1:d, l1 in l2:d# only the lower triangular part
    C[l1, l2] = sum(ev .* ev ./ (λ[l1] * ev + oneT) ./ (λ[l2] * ev + oneT))
  end
  LinAlg.copytri!(C, 'L') # copy to upper triangular part
  A_mul_Bt!(sub(H, 1:d^2, 1:d^2), scale(Φ2, vec(C)), Φ2)
  # (2, 1) block
  for l2 in 1:d, l1 in l2:d
    C[l1, l2] = sum(ev ./ (λ[l1] * ev + oneT) ./ (λ[l2] * ev + oneT))
  end
  LinAlg.copytri!(C, 'L') # copy to upper triangular part
  A_mul_Bt!(sub(H, d^2+1:2d^2, 1:d^2), scale(Φ2, vec(C)), Φ2)
  # d-by-d (2, 2) block
  for l2 in 1:d, l1 in l2:d
    C[l1, l2] = sum(oneT ./ (λ[l1] * ev + oneT) ./ (λ[l2] * ev + oneT))
  end
  LinAlg.copytri!(C, 'L') # copy to upper triangular part
  A_mul_Bt!(sub(H, d^2+1:2d^2, d^2+1:2d^2), scale(Φ2, vec(C)), Φ2)
  # make sure it's Hessian of the *negative* log-likehood
  LinAlg.copytri!(H, 'L') # copy to upper triangular part
  scale!(H, convert(T, 0.5))
end # function reml_fisher

function reml_fisher!{T <: AbstractFloat}(
  Σ::Union{Vector{Matrix{T}}, Tuple{Matrix{T}, Matrix{T}}},
  ev::Vector{T}
  )

  d = size(Σ[1], 1)
  H = zeros(T, 2d^2, 2d^2)
  reml_fisher!(H, Σ, ev)
  H
end

"""

    reml_eig(Y, V)

Extract eigen-decomposition of `V = (V[1], V[2])`.

# Input
- `Y`: `n x d` response variables.
- `V = (V[1], V[2])`: two `n x n` covariance matrices.

# Output
- `Yrot`: `n x d` rotated response variables `U' * Y`.
- `deval`: eigenvalues from eigen-decomposition `(deval, U) = eig(V[1], V[2])`.
- `loglconst`: constant `n*d*log(2π) + d*logdet(V2)` in 2.0log-likelihood.
"""
function reml_eig{T <: AbstractFloat}(
  Y::AbstractVecOrMat{T},
  V::Union{Vector{Matrix{T}}, Tuple{Matrix{T}, Matrix{T}}}
  )

  n, d = size(Y, 1), size(Y, 2)
  nd = n * d
  # (generalized)-eigendecomposition of (V1, V2)
  if isa(V[2], UniformScaling)
    deval, U = eig(V[1])
    loglconst = - nd * log(2.0π)
  elseif isdiag(V[2])
    deval, U = eig(V[1])
    loglconst = - nd * log(2.0π) - d * sum(log, diag(V[2]))
  else
    deval, U = eig(V[1], V[2])
    #scale!(U, 1.0 ./ sqrt(diag(U' * V[2] * U)))
    loglconst = - nd * log(2.0π) - d * logdet(V[2])
  end
  # corect negative eigenvalues due to roundoff error
  zeroT = convert(T, 0.0)
  map!(x -> max(x, zeroT), deval) # correct negative eigenvalues due to roundoff
  # rotate responses
  Yrot = At_mul_B(U, Y)
  # output
  Yrot, deval, loglconst
end # function reml_eig

# Set up MathProgBase interface

type TwoVarComp{T <: AbstractFloat} <: MathProgBase.AbstractNLPEvaluator
  Yrot::AbstractVecOrMat{T}
  ev::Vector{T}
  loglconst::T
  L::Vector{Matrix{T}}
  Σ::Vector{Matrix{T}}
  ∇Σ::Vector{T} # graident wrt Σs
  HΣ::Matrix{T} # Hessian wrt Σs
  HL::Matrix{T} # Hessian wrt Ls
end

function MathProgBase.initialize(dd::TwoVarComp,
  requested_features::Vector{Symbol})
  for feat in requested_features
    if !(feat in [:Grad, :Jac, :Hess])
      error("Unsupported feature $feat")
    end
  end
end # function MathProgBase.initialize

MathProgBase.features_available(dd::TwoVarComp) = [:Grad, :Jac, :Hess]
MathProgBase.eval_g(dd::TwoVarComp, g, x) = nothing
MathProgBase.jac_structure(dd::TwoVarComp) = Int[], Int[]
MathProgBase.eval_jac_g(dd::TwoVarComp, J, x) = nothing

function MathProgBase.eval_f{T}(dd::TwoVarComp, x::Vector{T})
  d = size(dd.L[1], 1)
  nparam = d * (d + 1)
  nparamhalf = div(nparam, 2)
  dd.L[1][trilind(dd.L[1])] = x[1:nparamhalf]
  dd.L[2][trilind(dd.L[2])] = x[nparamhalf+1:end]
  A_mul_Bt!(dd.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(dd.Σ[2], dd.L[2], dd.L[2])
  reml_objval(dd.Σ, dd.Yrot, dd.ev, dd.loglconst)
end # function MathProgBase.eval_f

function MathProgBase.eval_grad_f{T}(dd::TwoVarComp,
  grad_f::Vector{T}, x::Vector{T})
  d = size(dd.L[1], 1)
  nparam = d * (d + 1)
  nparamhalf = div(nparam, 2)
  dd.L[1][trilind(dd.L[1])] = x[1:nparamhalf]
  dd.L[2][trilind(dd.L[2])] = x[nparamhalf+1:end]
  A_mul_Bt!(dd.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(dd.Σ[2], dd.L[2], dd.L[2])
  reml_grad!(dd.∇Σ, dd.Σ, dd.Yrot, dd.ev)
  # chain rule for gradient wrt Cholesky factor
  chol_gradient!(sub(grad_f, 1:nparamhalf),
    dd.∇Σ[1:d^2], dd.L[1])
  chol_gradient!(sub(grad_f, nparamhalf+1:nparam),
    dd.∇Σ[d^2+1:end], dd.L[2])
end # function MathProgBase.eval_grad_f{

function MathProgBase.hesslag_structure(dd::TwoVarComp)
  d = size(dd.L[1], 1)
  nparam = d * (d + 1)
  ind2sub((nparam, nparam), trilind(nparam))
end # function MathProgBase.hesslag_structure

function MathProgBase.eval_hesslag{T}(dd::TwoVarComp, H::Vector{T},
  x::Vector{T}, σ::T, μ::Vector{T})
  d = size(dd.L[1], 1)
  nparam = d * (d + 1)
  nparamhalf = div(nparam, 2)
  dd.L[1][trilind(dd.L[1])] = x[1:nparamhalf]
  dd.L[2][trilind(dd.L[2])] = x[nparamhalf+1:end]
  A_mul_Bt!(dd.Σ[1], dd.L[1], dd.L[1])
  A_mul_Bt!(dd.Σ[2], dd.L[2], dd.L[2])
  reml_fisher!(dd.HΣ, dd.Σ, dd.ev)
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

"""
    reml_fs(Yrot, ev, loglconst; Σ0, maxiter, solver)

Fit variance component model using the Fisher's scoring algorithm. Data `vec(Y)`
is assumed to be normal with mean zero and covariance `Σ[1]⊗V[1] + Σ[2]⊗V[2]`.

# Input
- `Yrot`: rotated responses `U' * Y`, where `(ev, U) = eig(V[1], V[2])`.
- `ev`: eigenvalues from `(ev, U) = eig(V[1], V[2])`.
- `loglconst`: constant `n*d*log(2π) + d*logdet(V2)` in 2.0log-likelihood.

# Keyword argument
- `Σ0=(Σ0[1], Σ0[2])`: starting value for variance component parameters.
- `maxiter`: maximum number of iterations for nonlinear programming solver.
- `solver`: nonlinear programming solver. Default is `:Ipopt`; other supported
  solvers are ':Knitro'.
- `verbose`: logging information.

# Output
- `logl`: log-likelihood at `Σ = (Σ[1], Σ[2])`.
- `Σ=(Σ[1], Σ[2])`: variance component estimates.
- `Σse=(Σse[1], Σse[2])`: standard errors of variance component estimates.
- `Σcov`: `2d^2 x 2d^2` covariance matrix of variance component estimates.
"""
function reml_fs{T <: AbstractFloat}(
  Yrot::AbstractVecOrMat{T},
  ev::Vector{T},
  loglconst::T;
  Σ0::Union{Vector{Matrix{T}}, Tuple{Matrix{T}, Matrix{T}}} =
    [eye(T, size(Yrot, 2)) for i = 1:2],
  maxiter::Integer = 1000,
  solver::Symbol = :Ipopt,
  verbose::Bool = true
  )

  n, d = size(Yrot, 1), size(Yrot, 2)
  nd = n * d
  # number of optimization parameters in Cholesky factors
  nparam = d * (d + 1)
  nparamhalf = div(nparam, 2)

  # pre-allocate variables for optimization
  zeroT = convert(T, 0.0)
  L = Array(Matrix{T}, 2)
  L[1] = zeros(T, d, d)
  L[2] = zeros(T, d, d)
  Σ = Array(Matrix{T}, 2)
  Σ[1] = zeros(T, d, d)
  Σ[2] = zeros(T, d, d)
  Ltrilind = trilind(d, d)
  ∇Σ = zeros(T, 2d^2) # graident wrt Σs
  HΣ = zeros(T, 2d^2, 2d^2) # Hessian wrt Σs
  HL = zeros(T, nparam, nparam) # Hessian wrt Ls
  # data for the optimization problem
  dd = TwoVarComp(Yrot, ev, loglconst, L, Σ, ∇Σ, HΣ, HL)

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
  x0 = [vech(chol(Σ0[1], Val{:L}).data); vech(chol(Σ0[2], Val{:L}).data)]
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
  A_mul_Bt!(Σ[1], L[1], L[1])
  A_mul_Bt!(Σ[2], L[2], L[2])

  # standard errors
  Σcov = zeros(T, 2d^2, 2d^2)
  reml_fisher!(Σcov, Σ, ev)
  Σcov = inv(Σcov)
  Σse = deepcopy(Σ)
  copy!(Σse[1], sqrt(diag(sub(Σcov, 1:d^2, 1:d^2))))
  copy!(Σse[2], sqrt(diag(sub(Σcov, d^2+1:2d^2, d^2+1:2d^2))))

  # output
  maxlogl, Σ, Σse, Σcov
end # function reml_fs

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
function reml_mm{T <: AbstractFloat}(
  Yrot::AbstractVecOrMat{T},
  ev::Vector{T},
  loglconst::T;
  Σ0::Union{Vector{Matrix{T}}, Tuple{Matrix{T}, Matrix{T}}} =
    [eye(T, size(Yrot, 2)) for i = 1:2],
  maxiter::Integer = 10000,
  funtol::T = convert(T, 1e-8),
  verbose::Bool = true)

  # initialize algorithm
  # n = no. observations, d = no. categories
  n, d = size(Yrot, 1), size(Yrot, 2)
  nd = n * d
  Σ = deepcopy(Σ0)
  zeroT, oneT, halfT = zero(T), one(T), convert(T, 0.5)
  # update generalized eigen-decomposition
  λ, Φ = eig(Σ[1], Σ[2])
  # normalize columns of Φ to be orthonormal with respect to Σ[2]
  # scale!(Φ, 1.0 ./ sqrt(diag(Φ' * Σ[2] * Φ)))
  Wt = oneT ./ sqrt(ev * λ' + oneT)
  res = (Yrot * Φ) .* Wt
  Wt = oneT ./ sqrt(ev * λ' + oneT)
  res = (Yrot * Φ) .* Wt
  logl = halfT * (loglconst - n * logdet(Σ[2]) - sumabs2(res)) + sum(log, Wt)
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
      @inbounds dg[j] = sqrt(sum(ev ./ (λ[j] * ev + oneT)))
    end
    Whalf = res .* Wt
    scale!(sqrt(ev), Whalf)
    scale!(Whalf, λ .* dg)
    #W = sqrtm(Whalf' * Whalf) # produces imaginery eigenvalues due to precision
    #dg = 1.0 ./ dg
    #Σ[1] = scale(dg, inv(Φ))' * W * scale(dg, inv(Φ))
    # this approach is more numerical stable
    Whalfsvd = svdfact(Whalf)
    Σ[1] = scale(sqrt(Whalfsvd[:S]), Whalfsvd[:Vt]) * scale(oneT ./ dg, inv(Φ))
    Σ[1] = Σ[1]' * Σ[1]
    # update Σ2
    @inbounds for j = 1:d
      dg[j] = sqrt(sum(oneT ./ (λ[j] * ev + oneT)))
    end
    Whalf = res .* Wt
    scale!(Whalf, dg)
    # W = sqrtm(Whalf' * Whalf)
    # dg = 1.0 ./ dg
    # Σ[2] = scale(dg, inv(Φ))' * W * scale(dg, inv(Φ))
    # this approach is more numerical stable
    Whalfsvd = svdfact(Whalf)
    Σ[2] = scale(sqrt(Whalfsvd[:S]), Whalfsvd[:Vt]) * scale(oneT ./ dg, inv(Φ))
    Σ[2] = Σ[2]' * Σ[2]

    # update generalized eigen-decomposition
    λ, Φ = eig(Σ[1], Σ[2])
    # normalize columns of Φ to be orthonormal with respect to Σ[2]
    # scale!(Φ, 1.0 ./ sqrt(diag(Φ' * Σ[2] * Φ)))
    Wt = oneT ./ sqrt(ev * λ' + oneT)
    res = (Yrot * Φ) .* Wt

    # check convergence
    loglold = logl
    logl = halfT * (loglconst - n * logdet(Σ[2]) - sumabs2(res)) + sum(log, Wt)
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
  reml_fisher!(Σcov, Σ, ev)
  Σcov = inv(Σcov)
  Σse = deepcopy(Σ)
  copy!(Σse[1], sqrt(diag(sub(Σcov, 1:d^2, 1:d^2))))
  copy!(Σse[2], sqrt(diag(sub(Σcov, d^2+1:2d^2, d^2+1:2d^2))))

  # output
  logl, Σ, Σse, Σcov
end # function reml_mm

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
function heritability{T <: AbstractFloat}(
  Σ::Union{Vector{Matrix{T}}, Tuple{Matrix{T}}},
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

  map(x -> logpdf(vcm, x), vcobs)
end

#---------------------------------------------------------------------------#
# Evaluate gradient
#---------------------------------------------------------------------------#

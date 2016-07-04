module VarianceComponentTypeTest

using VarianceComponentModels, MathProgBase, Ipopt
using BaseTestNext

srand(123)

# generate data from a d-variate response variane component model
n = 100   # no. observations
d = 2     # no. categories
m = 2     # no. variance components
p = 2     # no. covariates
X = randn(n, p)
B = ones(p, d)
Σ = ntuple(x -> zeros(d, d), m)
for i in 1:m
  Σi = randn(d, d)
  copy!(Σ[i], Σi' * Σi)
end
## make the first variance component 0 matrix
#fill!(Σ[1], 0.0)
V = ntuple(x -> zeros(n, n), m)
for i = 1:m-1
  Vi = randn(n, 50)
  copy!(V[i], Vi * Vi')
end
copy!(V[m], eye(n))
# form Ω
Ω = zeros(n*d, n*d)
for i = 1:m
  Ω += kron(Σ[i], V[i])
end
Ωchol = cholfact(Ω)
Y = X * B + reshape(Ωchol[:L] * randn(n*d), n, d)

info("Forming variance component model and data")
vcdata = VarianceComponentVariate(Y, X, V)
vcmodel = VarianceComponentModel(vcdata)

info("Pre-compute (generalized) eigen-decomposition and rotate data")
vcdatarot = TwoVarCompVariateRotate(vcdata)
vcmodelrot = TwoVarCompModelRotate(vcmodel)

info("Evaluate log-pdf")
#@code_warntype logpdf(vcmodelrot, vcdatarot)
@inferred logpdf(vcmodelrot, vcdatarot)
@test logpdf(vcmodel, vcdata) == logpdf(vcmodelrot, vcdatarot)
# @test (logpdf(vcmodelrot, [vcdatarot vcdatarot; vcdatarot vcdatarot]) -
#   logpdf(vcmodel, [vcdata vcdata; vcdata vcdata])) ≈ 0.0

info("Evaluate gradient")
∇ = zeros(nmeanparams(vcmodel) + 2d^2)
#@code_warntype gradient!(∇, vcmodelrot, vcdatarot)
@inferred gradient!(∇, vcmodelrot, vcdatarot)
@test vecnorm(gradient(vcmodel, vcdata) - gradient(vcmodelrot, vcdatarot)) ≈ 0.0

info("Evaluate Fisher information matrix")
H = zeros(nmeanparams(vcmodel) + 2d^2, nmeanparams(vcmodel) + 2d^2)
#@code_warntype fisher!(H, vcmodelrot, vcdatarot)
@inferred fisher!(H, vcmodelrot, vcdatarot)
@test vecnorm(fisher(vcmodel, vcdata) - fisher(vcmodelrot, vcdatarot)) ≈ 0.0

info("Update mean parameters from variance component parameters (no constraints)")
vcmmean = deepcopy(vcmodel)
#@code_warntype update_meanparam!(vcmmean, vcdatarot)
#@inferred update_meanparam!(vcmodel, vcdatarot)
update_meanparam!(vcmmean, vcdatarot)
@show vcmmean.B

info("Update mean parameters from variance component parameters (linear equality constraints)")
vcmmean = deepcopy(vcmodel)
#@code_warntype update_meanparam!(vcmmean, vcdatarot, [1.0 -1.0 0 0], '=', 0.0)
update_meanparam!(vcmmean, vcdatarot, [1.0 -1.0 zeros(1, p*d-2)], '=', 0.0)
@show vcmmean.B
@test vcmmean.B[1] ≈ vcmmean.B[2]

info("Update mean parameters from variance component parameters (box constraints)")
vcmmean = deepcopy(vcmodel)
# @code_warntype update_meanparam!(vcmmean, vcdatarot, zeros(0, p * d), Vector{Char}(0),
#   Vector{Float64}(0), 0.0, 1.0, IpoptSolver(print_level = 0), zeros(n * d, p * d), zeros(n * d),
#   zeros(p * d, p * d), zeros(p * d))
update_meanparam!(vcmmean, vcdatarot, Matrix{Float64}(0, p * d), Vector{Char}(0),
  Vector{Float64}(0), 0.0, 1.0)
@show vcmmean.B
@test all(vcmmean.B .≥ 0.0)
@test all(vcmmean.B .≤ 1.0)

info("Find MLE using Fisher scoring")
vcmfs = deepcopy(vcmodel)
logl_fs, _, _, Σcov_fs, Bse_fs, = mle_fs!(vcmfs, vcdatarot; solver = :Ipopt)
@show vcmfs.B, Bse_fs, B

info("Find MLE using MM algorithm")
vcmm = deepcopy(vcmodel)
# @code_warntype mle_mm!(vcmm, vcdatarot;
#   maxiter = 10000,
#   funtol = 1e-8,
#   verbose = true,
#   A = zeros(0, p * d),
#   sense = '=',
#   b = 0.0,
#   lb = -Inf,
#   ub = Inf,
#   qpsolver = IpoptSolver(print_level = 0)
# )
#@inferred mle_mm!(vcmm, vcdatarot)
logl_mm, _, _, Σcov_mm = mle_mm!(vcmm, vcdatarot)
@test abs(logl_fs - logl_mm) / (abs(logl_fs) + 1.0) < 1.0e-4

info("Find MLE using Fisher scoring (linear equality + box constraints)")
vcmfs = deepcopy(vcmodel)
logl_fs, _, _, Σcov_fs, Bse_fs, = mle_fs!(vcmfs, vcdatarot; solver = :Ipopt,
  A = [1.0 -1.0 zeros(1, p*d-2)], sense = '=', b = 0.0, lb = 0.0, ub = 1.0)
@show vcmfs.B
@test vcmfs.B[1] ≈ vcmfs.B[2]
@test all(vcmfs.B .≥ 0.0)
@test all(vcmfs.B .≤ 1.0)

info("Find MLE using MM algorithm (linear equality + box constraints)")
vcmm = deepcopy(vcmodel)
logl_mm, _, _, Σcov_mm = mle_mm!(vcmm, vcdatarot; A = [1.0 -1.0 zeros(1, p*d-2)],
  sense = '=', b = 0.0, lb = 0.0, ub = 1.0)
@show vcmm.B
@test vcmm.B[1] ≈ vcmm.B[2]
@test all(vcmm.B .≥ 0.0)
@test all(vcmm.B .≤ 1.0)
@test abs(logl_fs - logl_mm) / (abs(logl_fs) + 1.0) < 1.0e-4

info("Heritability estimation")
h, h_se = heritability(vcmfs.Σ, Σcov_fs)
@show h, h_se
@test all(0.0 .≤ h .≤ 1.0)

info("test fit_mle")
vcmmle = deepcopy(vcmodel)
logl_mle, _, _, Σcov_mle, Bse_mle, = fit_mle!(vcmmle, vcdata; algo = :FS)
@show vcmfs.B, Bse_fs, B

info("test fit_reml")
vcmreml = deepcopy(vcmodel)
logl_reml, _, _, Σcov_reml, Bse_reml, = fit_reml!(vcmreml, vcdata; algo = :FS)
@show vcmreml.B, Bse_reml, B

end # module VarianceComponentTypeTest

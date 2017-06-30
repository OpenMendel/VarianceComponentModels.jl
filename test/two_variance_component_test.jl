module TwoVarianceComponentTest

using VarianceComponentModels, MathProgBase, Ipopt, Base.Test

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

info("Forming VarianceComponentModel from data")
@inferred VarianceComponentVariate(Y, X, V)
vcdata  = VarianceComponentVariate(Y, X, V)
vcmodel = VarianceComponentModel(vcdata)

info("Pre-compute eigen-decomposition and rotate data")
vcdatarot  = TwoVarCompVariateRotate(vcdata)
vcmodelrot = TwoVarCompModelRotate(vcmodel)

info("Evaluate log-pdf")
#@code_warntype logpdf(vcmodelrot, vcdatarot)
@inferred logpdf(vcmodelrot, vcdatarot)
@test logpdf(vcmodel, vcdata) == logpdf(vcmodelrot, vcdatarot)
@test (logpdf(vcmodelrot, [vcdatarot vcdatarot; vcdatarot vcdatarot]) -
  logpdf(vcmodel, [vcdata vcdata; vcdata vcdata])) ≈ 0.0

info("Evaluate gradient")
∇ = zeros(2d^2)
#@code_warntype gradient!(∇, vcmodelrot, vcdatarot)
@inferred gradient!(∇, vcmodelrot, vcdatarot)
@test vecnorm(gradient(vcmodel, vcdata) - gradient(vcmodelrot, vcdatarot)) ≈ 0.0
@test vecnorm(gradient(vcmodel, vcdata) - gradient(vcmodel, vcdatarot)) ≈ 0.0
@test vecnorm(gradient(vcmodel, [vcdata vcdata]) -
  2.0gradient(vcmodel, vcdata)) ≈ 0.0
@test vecnorm(gradient(vcmodel, [vcdata vcdata]) -
  gradient(vcmodelrot, [vcdatarot vcdatarot])) ≈ 0.0

info("Evaluate Fisher information matrix of Σ")
H = zeros(2d^2, 2d^2)
#@code_warntype fisher!(H, vcmodelrot, vcdatarot)
@inferred fisher_Σ!(H, vcmodelrot, vcdatarot)
@test vecnorm(fisher_Σ(vcmodel, vcdata) - fisher_Σ(vcmodelrot, vcdatarot)) ≈ 0.0
@test vecnorm(fisher_Σ(vcmodel, vcdata) - fisher_Σ(vcmodel, vcdatarot)) ≈ 0.0
@test vecnorm(fisher_Σ(vcmodel, [vcdata vcdata]) -
  2fisher_Σ(vcmodel, vcdata)) ≈ 0.0
@test vecnorm(fisher_Σ(vcmodel, [vcdata vcdata]) -
  fisher_Σ(vcmodelrot, [vcdatarot vcdatarot])) ≈ 0.0


info("Evaluate Fisher information matrix of B")
H = zeros(p * d, p * d)
#@code_warntype fisher_B!(H, vcmodelrot, vcdatarot)
@inferred fisher_B!(H, vcmodelrot, vcdatarot)
@test vecnorm(fisher_B(vcmodel, vcdata) - fisher_B(vcmodelrot, vcdatarot)) ≈ 0.0
@test vecnorm(fisher_B(vcmodel, vcdata) - fisher_B(vcmodel, vcdatarot)) ≈ 0.0
@test vecnorm(fisher_B(vcmodel, [vcdata vcdata]) -
  2.0fisher_B(vcmodel, vcdata)) ≈ 0.0
@test vecnorm(fisher_B(vcmodel, [vcdata vcdata]) -
  fisher_B(vcmodelrot, [vcdatarot vcdatarot])) ≈ 0.0

info("Find MLE using Fisher scoring")
vcmfs = deepcopy(vcmodel)
#@code_warntype mle_fs!(vcmfs, vcdatarot; solver = :Ipopt)
#@inferred mle_fs!(vcmfs, vcdatarot; solver = :Ipopt)
logl_fs, _, _, Σcov_fs, Bse_fs, = mle_fs!(vcmfs, vcdatarot; solver = :Ipopt)
logl_fs_array, = mle_fs!(vcmfs, [vcdatarot vcdatarot]; solver = :Ipopt)
@show vcmfs.B
@show Bse_fs
@show B

info("Find MLE using MM algorithm")
vcmmm = deepcopy(vcmodel)
#@code_warntype mle_mm!(vcmm, vcdatarot)
#@inferred mle_mm!(vcmmm, vcdatarot)
logl_mm, _, _, _, Bse_mm, = mle_mm!(vcmmm, vcdatarot)
@test abs(logl_fs - logl_mm) / (abs(logl_fs) + 1.0) < 1.0e-4
@show vcmmm.B
@show Bse_mm
@show B
vcmmm = deepcopy(vcmodel)
logl_mm_array, = mle_mm!(vcmmm, [vcdatarot vcdatarot])
@test abs(logl_fs_array - logl_mm_array) / (abs(logl_fs_array) + 1.0) < 1.0e-4

info("Find MLE using Fisher scoring (linear equality + box constraints)")
vcmfs       = deepcopy(vcmodel)
vcmfs.A     = [1.0 -1.0 zeros(1, p*d-2)]
vcmfs.sense = '='
vcmfs.b     = 0.0
vcmfs.lb    = 0.0
vcmfs.ub    = 1.0
logl_fs, _, _, Σcov_fs, Bse_fs, = mle_fs!(vcmfs, vcdatarot; solver = :Ipopt, qpsolver = :Ipopt)
@show vcmfs.B
@test vcmfs.B[1] ≈ vcmfs.B[2]
@test all(vcmfs.B .≥ 0.0)
@test all(vcmfs.B .≤ 1.0)

info("Find MLE using MM algorithm (linear equality + box constraints)")
vcmm = deepcopy(vcmodel)
vcmm.A = [1.0 -1.0 zeros(1, p*d-2)]
vcmm.sense = '='
vcmm.b = 0.0
vcmm.lb = 0.0
vcmm.ub = 1.0
logl_mm, _, _, Σcov_mm = mle_mm!(vcmm, vcdatarot; qpsolver = :Ipopt)
@show vcmm.B
@test vcmm.B[1] ≈ vcmm.B[2]
@test all(vcmm.B .≥ 0.0)
@test all(vcmm.B .≤ 1.0)
@test abs(logl_fs - logl_mm) / (abs(logl_fs) + 1.0) < 1.0e-4

info("Heritability estimation")
h, h_se = heritability(vcmfs.Σ, Σcov_fs)
@show h, h_se
@test all(0.0 .≤ h .≤ 1.0)

info("test fit_mle (FS)")
vcmmle = deepcopy(vcmodel)
logl_mle, _, _, Σcov_mle, Bse_mle, = fit_mle!(vcmmle, vcdata; algo = :FS)
@show vcmmle.B, Bse_mle, B
vcmle = deepcopy(vcmodel)
fit_mle!(vcmmle, [vcdata vcdata]; algo = :FS)

info("test fit_mle (MM)")
vcmmle = deepcopy(vcmodel)
logl_mle, _, _, Σcov_mle, Bse_mle, = fit_mle!(vcmmle, vcdata; algo = :MM)
@show vcmmle.B, Bse_mle, B

info("test fit_reml (FS)")
vcmreml = deepcopy(vcmodel)
logl_reml, _, _, Σcov_reml, Bse_reml, = fit_reml!(vcmreml, vcdata; algo = :FS)
@show vcmreml.B, Bse_reml, B

info("test fit_reml (MM)")
vcmreml = deepcopy(vcmodel)
logl_reml, _, _, Σcov_reml, Bse_reml, = fit_reml!(vcmreml, vcdata; algo = :MM)
@show vcmreml.B, Bse_reml, B

end # module VarianceComponentTypeTest

module VarianceComponentTypeTest

using VarianceComponentModels
using BaseTestNext

srand(123)

# generate data from a d-variate response variane component model
n = 100   # no. observations
d = 3     # no. categories
m = 2     # no. variance components
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
Y = reshape(Ωchol[:L] * randn(n*d), n, d)

info("Forming variance component model and data")
vcdata = VarianceComponentVariate(Y, Float64[], V)
vcmodel = VarianceComponentModel(Float64[], Σ)

info("Pre-compute (generalized) eigen-decomposition and rotate data")
vcdatarot = TwoVarCompVariateRotate(vcdata)
vcmodelrot = TwoVarCompModelRotate(vcmodel)

info("Evaluate log-pdf")
@test logpdf(vcmodel, vcdata) == logpdf(vcmodelrot, vcdatarot)
@test vecnorm(logpdf(vcmodel, [vcdata vcdata; vcdata vcdata]) -
  logpdf(vcmodelrot, [vcdatarot vcdatarot; vcdatarot vcdatarot])) < 1.0e-8

info("Evaluate gradient")
@test vecnorm(gradient(vcmodel, vcdata) - gradient(vcmodelrot, vcdatarot)) ≈ 0.0

info("Evaluate Fisher information matrix")
@test vecnorm(fisher(vcmodel, vcdata) - fisher(vcmodelrot, vcdatarot)) ≈ 0.0

info("Find MLE using Fisher scoring")
vcmfs = deepcopy(vcmodel)
logl_fs, _, Σcov_fs = mle_fs!(vcmfs, vcdatarot; solver = :Ipopt)

info("Find MLE using MM algorithm")
vcmmm = deepcopy(vcmodel)
logl_mm, _, Σcov_mm = mle_mm!(vcmmm, vcdatarot)
@test abs(logl_fs - logl_mm) / (abs(logl_fs) + 1.0) < 1.0e-4

info("Heritability estimation")
h, h_se = heritability(vcmfs.Σ, Σcov_fs)
@show h, h_se
@test all(0.0 .≤ h .≤ 1.0)

end # module VarianceComponentTypeTest

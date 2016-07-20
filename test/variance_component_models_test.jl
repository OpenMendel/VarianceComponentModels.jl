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

info("VarianceComponentModel type")
#@code_warntype VarianceComponentModel(B, Σ)
@inferred VarianceComponentModel(B, Σ)
#@code_warntype VarianceComponentModel(Σ)
@inferred VarianceComponentModel(Σ)
vcmodel = VarianceComponentModel(Σ)
@test isempty(vcmodel.B)

info("TwoVarCompModelRotate type")
F = eigfact(Σ[1], Σ[2])
#@code_warntype TwoVarCompModelRotate(B * F.vectors, F.values, F.vectors, logdet(Σ[2]))
@inferred TwoVarCompModelRotate(B * F.vectors, F.values, F.vectors, logdet(Σ[2]))

info("VarianceComponentVariate type")
#@code_warntype VarianceComponentVariate(Y, V)
@inferred VarianceComponentVariate(Y, V)
vcdata = VarianceComponentVariate(Y, V)
@test isempty(vcdata.X)

info("Forming VarianceComponentModel from data")
vcdata = VarianceComponentVariate(Y, X, V)
#@code_warntype VarianceComponentVariate(Y, X, V)
@inferred VarianceComponentVariate(Y, X, V)
vcdata = VarianceComponentVariate(Y, X, V)
#@code_warntype VarianceComponentModel(vcdata)
@inferred VarianceComponentModel(vcdata)
vcmodel = VarianceComponentModel(vcdata)
@test vecnorm(vcmodel.B) ≈ 0.0
@test vecnorm(vcmodel.Σ[1] - eye(d)) ≈ 0.0
@test vecnorm(vcmodel.Σ[2] - eye(d)) ≈ 0.0

info("Pre-compute eigen-decomposition and rotate data")
#@code_warntype TwoVarCompVariateRotate(vcdata)
@inferred TwoVarCompVariateRotate(vcdata)
vcdatarot = TwoVarCompVariateRotate(vcdata)
#@code_warntype TwoVarCompModelRotate(vcmodel)
@inferred TwoVarCompModelRotate(vcmodel)
vcmodelrot = TwoVarCompModelRotate(vcmodel)

info("Test generalized eigen-decomposition")
vcdatacopy = deepcopy(vcdata)
scale!(vcdatacopy.V[2], 2.0) # V[2] = 2I
vcdatacopyrot = TwoVarCompVariateRotate(vcdatacopy)
@test vcdatacopyrot.logdetV2 ≈ n * log(2.0) # det(V[2]) = 2^n

info("VarianceComponentModel from TwoVarCompVariateRotate")
@inferred VarianceComponentModel(vcdatarot)
vcmodelfromrot = VarianceComponentModel(vcdatarot)
@test size(vcmodelfromrot.B, 1) == p
@test size(vcmodelfromrot.B, 2) == d
@test length(vcmodelfromrot.Σ) == m

info("VarianceComponentAuxData from VarianceComponentVariate")
vcobsaux1 = VarianceComponentAuxData(vcdata)
vcobsaux2 = VarianceComponentAuxData(vcdatarot)
@test size(vcobsaux1.res) == size(vcobsaux2.res)
@test size(vcobsaux1.Xwork) == size(vcobsaux2.Xwork)
@test size(vcobsaux1.ywork) == size(vcobsaux2.ywork)
@test size(vcobsaux1.obswt) == size(vcobsaux2.obswt)

info("Query functions")
@test eltype(vcmodel) == eltype(B)
@test eltype(vcdata) == eltype(Y)
@test eltype(vcmodelrot) == eltype(vcmodel)
@test eltype(vcdatarot) == eltype(vcdata)
@test length(vcmodel) == d
@test length(vcdata) == d
@test length(vcmodelrot) == d
@test length(vcdatarot) == d
@test size(vcdata) == (n, d)
@test size(vcdatarot) == (n, d)
@test nvarcomps(vcmodel) == m
@test nvarcomps(vcdata) == m
@test nvarcomps(vcmodelrot) == m
@test nvarcomps(vcdatarot) == m
@test nmeanparams(vcmodel) == p * d
@test nmeanparams(vcmodelrot) == p * d
@test nmeanparams(vcdata) == p * d
@test nmeanparams(vcdatarot) == p * d
@test nvarparams(vcmodel) == m * binomial(d + 1, 2)
@test nvarparams(vcmodelrot) == m * binomial(d + 1, 2)
@test nparams(vcmodel) == p * d + m * binomial(d + 1, 2)
@test nparams(vcmodelrot) == p * d + m * binomial(d + 1, 2)

info("Mean, covariance, and residual of model")
#@code_warntype mean(vcmodel, vcdata)
@inferred mean(vcmodel, vcdata)
μ = mean(vcmodel, vcdata)
@test size(μ, 1) == n
@test size(μ, 2) == d
#@code_warntype cov(vcmodel, vcdata)
@inferred cov(vcmodel, vcdata)
Ω = cov(vcmodel, vcdata)
@test size(Ω, 1) == n * d
@test size(Ω, 2) == n * d
#@code_warntype residual(vcmodel, vcdata)
@inferred residual(vcmodel, vcdata)
res = residual(vcmodel, vcdata)
@test size(res, 1) == n
@test size(res, 2) == d
@inferred residual(vcmodelrot, vcdatarot)
resrot = residual(vcmodelrot, vcdatarot)
@test size(resrot, 1) == n
@test size(resrot, 2) == d

end # module VarianceComponentTypeTest

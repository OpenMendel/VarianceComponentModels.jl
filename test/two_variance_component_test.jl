module TwoVarianceComponentTest

using VarianceComponentModels
using BaseTestNext

srand(123)

# generate data from a d-variate response variane component model
n = 100   # no. observations
d = 2     # no. categories
m = 2     # no. variance components
Σ = Array(Matrix{Float64}, m)
for i = 1:m
  Σi = randn(d, d)
  Σ[i] = Σi' * Σi
end
## make the first variance component 0 matrix
#fill!(Σ[1], 0.0)
V = Array(Matrix{Float64}, m)
for i = 1:m-1
  Vi = randn(n, 50)
  V[i] = Vi * Vi'
end
V[m] = eye(n)
# form Ω
Ω = zeros(n*d, n*d)
for i = 1:m
  Ω += kron(Σ[i], V[i])
end
Ωchol = cholfact(Ω)
Y = reshape(Ωchol[:L] * randn(n*d), n, d)

@testset "fit vc and estimate heritability" begin
info("Pre-compute (generalized) eigen-decomposition")
# pre-compute the (generalized) eigen-decomposition
Yrot, deval, loglconst = reml_eig(Y, V)
info("Fit variance component model using Fisher scoring algorithm")
# fit variance component model using Fisher scoring
logl_fs, Σhat_fs, Σse_fs, Σcov_fs = reml_fs(Yrot, deval, loglconst; solver = :Ipopt)
@test vecnorm(reml_grad(Σhat_fs, Yrot, deval)) < 1.0e-3
info("Fit variance component model using MM algorithm")
# fit variance component model using MM algorithm
logl_mm, Σhat_mm, Σse_mm, Σcov_mm = reml_mm(Yrot, deval, loglconst)
#@test vecnorm(reml_grad(Σhat_mm, Yrot, deval)) < 1.0e-2
@test abs(logl_fs - logl_mm) / (abs(logl_fs) + 1.0) < 1.0e-3
info("Heritability estimation")
# heritability estimation
h, h_se = heritability(Σhat_fs, Σcov_fs)
@show h, h_se
@test all(0.0 .≤ h .≤ 1.0)
end

end # module MultivariateCalculusTest

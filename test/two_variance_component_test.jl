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

Yrot, deval, loglconst = reml_eig(Y, V)
logl_fs, Σhat_fs = reml_fs(Yrot, deval, loglconst; solver = :Ipopt)
@test vecnorm(reml_grad(Σhat_fs, Yrot, deval)) < 1.0e-3
logl_fs, Σhat_mm = reml_mm(Yrot, deval, loglconst)
@test vecnorm(reml_grad(Σhat_mm, Yrot, deval)) < 1.0e-2

end # module MultivariateCalculusTest

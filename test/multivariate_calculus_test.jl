module MultivariateCalculusTest

using VarianceComponentModels
using BaseTestNext

srand(123)

# test vech
@testset "vech" begin
  a = reshape(1:9, 3, 3)
  vecha = vech(a)
  @test vecha[1] == a[1, 1]
  @test vecha[2] == a[2, 1]
  @test vecha[3] == a[3, 1]
  @test vecha[4] == a[2, 2]
  @test vecha[5] == a[3, 2]
  @test vecha[6] == a[3, 3]
  a = 1.0
  @test vech(a) == a
end

# test trilind
@testset "trilind" begin
  n = 3
  A = randn(n, n)
  @test vecnorm(vech(A) - A[trilind(A)]) ≈ 0.0
end

# test triuind
@testset "triuind" begin
  n = 2
  idx = triuind(n, n)
  @test length(idx) == binomial(n + 1, 2)
  @test idx[1] == 1
  @test idx[2] == 3
  @test idx[3] == 4
  idx = triuind(n)
  @test length(idx) == binomial(n + 1, 2)
  @test idx[1] == 1
  @test idx[2] == 3
  @test idx[3] == 4
  idx = triuind(randn(n, n))
  @test length(idx) == binomial(n + 1, 2)
  @test idx[1] == 1
  @test idx[2] == 3
  @test idx[3] == 4
  idx = triuind(randn(n, n), 1)
  @test length(idx) == binomial(n, 2)
  @test idx[1] == 3
end

# test commutation
@testset "commutation" begin
  m, n = 3, 2
  A = randn(m, n)
  @test vecnorm(commutation(m, n) * vec(A) - vec(A')) ≈ 0.0
  @test vecnorm(commutation(A) * vec(A) - vec(A')) ≈ 0.0
  @test vecnorm(spcommutation(m, n) * vec(A) - vec(A')) ≈ 0.0
  @test vecnorm(spcommutation(A) * vec(A) - vec(A')) ≈ 0.0
  B = randn(m, m)
  @test vecnorm(commutation(m) * vec(B) - vec(B')) ≈ 0.0
  @test vecnorm(commutation(eltype(B), m) * vec(B) - vec(B')) ≈ 0.0
  @test vecnorm(spcommutation(eltype(B), m) * vec(B) - vec(B')) ≈ 0.0
end

# test duplication
@testset "duplication" begin
  n = 3
  A = randn(n, n)
  A = 0.5(A + A') # symmetrize
  @test vecnorm(duplication(n) * vech(A) - vec(A)) ≈ 0.0
  @test vecnorm(duplication(A) * vech(A) - vec(A)) ≈ 0.0
  @test vecnorm(spduplication(n) * vech(A) - vec(A)) ≈ 0.0
  @test vecnorm(spduplication(A) * vech(A) - vec(A)) ≈ 0.0
end

# test chol_gradient
@testset "chol_gradient" begin
  n = 3
  A = randn(n, n)
  A = 0.5(A + A') # symmetrize
  L = tril(randn(n, n))
  # calculate gradient wrt L using chol_gradient
  dL1 = chol_gradient(vec(A), L)
  # alternative way to calculate gradient wrt L
  dL2 = 2.0vech((A * L) + L' * A' - diagm(diag(A * L)))
  @test vecnorm(dL1 - dL2) ≈ 0.0
end

# test kron_gradient
@testset "kron_gradient" begin
  n, q, p, r = 2, 3, 4, 5
  X = randn(n, q)
  Y = randn(p, r)
  M = kron(X, Y)
  # calculate gradient wrt X using kron_gradient
  dX1 = kron_gradient(vec(M), Y, n, q)
  # alternative way to calculate gradient wrt X
  dX2 = zeros(X)
  for j = 1:q
    for i = 1:n
      dX2[i, j] = vecdot(Y, M[(i-1)*p+1:i*p, (j-1)*r+1:j*r])
    end
  end
  @test_approx_eq_eps vecnorm(dX1 - vec(dX2)) 0.0 1.0e-8
end

# test duplication
@testset "kronaxpy" begin
  m, n, p, q = 3, 4, 5, 6
  A = randn(m, n)
  X = randn(p, q)
  Y = zeros(m * p, n * q)
  kronaxpy!(A, X, Y)
  @test vecnorm(Y - kron(A, X)) ≈ 0.0
end

# test bump_diagonal
@testset "bump_diag" begin
  A = randn(3, 3)
  B = copy(A)
  bump_diagonal!(B, 1e-3)
  @test B[1, 1] - A[1, 1] ≈ 1e-3
  @test B[2, 2] - A[2, 2] ≈ 1e-3
  @test B[3, 3] - A[3, 3] ≈ 1e-3
end

# test clamp_diagonal
@testset "clamp_diag" begin
  n = 3
  A = randn(n, n)
  B = copy(A)
  clamp_diagonal!(B, 0.0, 1.0)
  @test all(0.0 .≤ diag(B) .≤ 1.0)
end


end # module MultivariateCalculusTest

module MultivariateCalculusTest

using MultivariateCalculus
using Base.Test

# test vech
a = reshape(1:15, 5, 3)
@show a
@show vech(a)

# test commutation
@show commutation(3, 2)
@show spcommutation(3, 2)
@show duplication(3)
@show spduplication(3)

# test chol_gradient
d = 3
A = randn(d, d); A = 0.5(A + A')
L = tril(randn(d, d))
@show 2.0vech((A * L) + L' * A' - diagm(diag(A * L)))
@show chol_gradient(vec(A), L)

srand(123)
n = 4; d = 3; nd = n * d;
A = randn(nd, nd); A = 0.5(A + A');
B = randn(nd, nd); B = 0.5(B + B');
#copy!(B, A) # make B same as A
V1 = randn(n, n); V1 = 0.5(V1 + V1');
V2 = randn(n, n); V2 = 0.5(V2 + V2');
result1 =
  MultivariateCalculus.kron_gradient(
  MultivariateCalculus.kron_gradient(kron(A, B), V1, d, d)', V1, d, d)'
result2 =
  kron(speye(d^2), vec(V1)') *
  kron(speye(d), MultivariateCalculus.spcommutation(d, n), speye(n)) *
  kron(A, B) *
  kron(speye(d), MultivariateCalculus.spcommutation(n, d), speye(n)) *
  kron(speye(d^2), vec(V1))
@show result1
@show result2
# (2, 1) block
for k2 in 1:d, j2 in 1:d, k1 in 1:d, j1 in 1:d
  println(trace(V1 * A[(k1-1)*n+1:k1*n, (k2-1)*n+1:k2*n] * V1 * B[(j2-1)*n+1:j2*n, (j1-1)*n+1:j1*n]))
end

end # module MultivariateCalculusTest

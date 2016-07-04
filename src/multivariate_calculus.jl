# module MultivariateCalculus

import Base.LowerTriangular
export vech, trilind, triuind,
  commutation, spcommutation,
  duplication, spduplication,
  chol_gradient, chol_gradient!,
  kron_gradient, kron_gradient!,
  kronaxpy!

"""
Vectorize the lower triangular part of a matrix.
"""
function vech(a::AbstractVecOrMat)
  # a is a vector
  if ndims(a) == 1
    return similar(a)
  end
  # a is a matrix
  m, n = size(a)
  out = similar(a, convert(Int, (2m - n + 1) * n / 2))
  ooffset, aoffset = 1, 1
  for j = 1:n
    len = m - j + 1 # no. elements to copy in column j
    copy!(out, ooffset, a, aoffset, len)
    ooffset += m - j + 1
    aoffset += m + 1
  end
  out
end

"""
    commutation(type, m[, n])

Create the `mn x mn` commutation matrix `K`, defined by
`K * vec(A) = vec(A')` for any `m x n` matrix A.
"""
function commutation(t::Type, m::Integer, n::Integer)
  ((m < 0) || (n < 0)) && throw(ArgumentError("invalid Array dimensions"))
  mn = m * n
  reshape(kron(vec(eye(t, m)), eye(t, n)), mn, mn)
end
commutation(m::Integer, n::Integer) = commutation(Float64, m, n)
commutation(t::Type, m::Integer) = commutation(t, m, m)
commutation(m::Integer) = commutation(m, m)
commutation(M::AbstractMatrix) = commutation(eltype(M), size(M, 1), size(M, 2))

"""
    spcommutation(type, m[, n])

Create the sparse `mn`-by-`mn` commutation matrix `K`, defined by
`K * vec(A) = vec(A')` for any `m x n` matrix A.
"""
function spcommutation(t::Type, m::Integer, n::Integer)
  ((m < 0) || (n < 0)) && throw(ArgumentError("invalid Array dimensions"))
  mn = m * n
  reshape(kron(vec(speye(t, m)), speye(t, n)), mn, mn)
end
spcommutation(m::Integer, n::Integer) = spcommutation(Float64, m, n)
spcommutation(t::Type, m::Integer) = spcommutation(t, m, m)
spcommutation(m::Integer) = spcommutation(m, m)
spcommutation(M::AbstractMatrix) = spcommutation(eltype(M), size(M, 1), size(M, 2))

"""
    trilind(m, n,[ k])

Linear indices of the lower triangular part of an `m x n` array.
"""
function trilind(m::Integer, n::Integer, k::Integer)
  find(tril(trues(m, n), k))
end
function trilind(m::Integer, n::Integer)
  find(tril(trues(m, n)))
end
function trilind(m::Integer)
  find(tril(trues(m, m)))
end
trilind(M::AbstractArray) = trilind(size(M, 1), size(M, 2))
trilind(M::AbstractArray, k::Integer) = trilind(size(M, 1), size(M, 2), k)

"""
    triuind(m, n,[ k])

Linear indices of the upper triangular part of an `m x n` array.
"""
function triuind(m::Integer, n::Integer, k::Integer)
  find(triu(trues(m, n), k))
end
function triuind(m::Integer, n::Integer)
  find(triu(trues(m, n)))
end
function triuind(m::Integer)
  find(triu(trues(m, m)))
end
triuind(M::AbstractArray) = triuind(size(M, 1), size(M, 2))
triuind(M::AbstractArray, k::Integer) = triuind(size(M, 1), size(M, 2), k)

"""
    spduplication(type, n)

Create the sparse `n^2 x n(n+1)/2` duplication matrix, defined by
`D * vech(A) = vec(A)` for any symmetric matrix.
"""
function spduplication(t::Type, n::Integer)
  imatrix = zeros(typeof(n), n, n)
  imatrix[trilind(n, n)] = 1:binomial(n + 1, 2)
  imatrix = imatrix + tril(imatrix, -1).'
  sparse(1:n^2, vec(imatrix), one(t))
end
spduplication(n::Integer) = spduplication(Float64, n)
spduplication(M::AbstractMatrix) = spduplication(eltype(M), size(M, 1))
duplication(t::Type, n::Integer) = full(spduplication(t, n))
duplication(n::Integer) = duplication(Float64, n)
duplication(M::AbstractMatrix) = duplication(eltype(M), size(M, 1))

"""
    kron_gradient!(g, dM, Y)

Compute the gradient `d / d vec(X)` from a vector of derivatives `dM` where
`M=X⊗Y`, `n, q = size(X)`, and `p, r = size(Y)`.
"""
function kron_gradient!{T <: Real}(g::VecOrMat{T}, dM::VecOrMat{T},
  Y::Matrix{T}, n::Integer, q::Integer)
  p, r = size(Y)
  A_mul_B!(g, kron(speye(n * q), vec(Y)'),
    (kron(speye(q), spcommutation(n, r), speye(p)) * dM))
end
function kron_gradient{T <: Real}(dM::VecOrMat{T}, Y::Matrix{T},
  n::Integer, q::Integer)
  if ndims(dM) == 1
    g = zeros(T, n * q)
  else
    g = zeros(T, n * q, size(dM, 2))
  end
  kron_gradient!(g, dM, Y, n, q)
end

"""
    chol_gradient!(g, dM, L)

Compute the gradient `d / d vech(L)` from a vector of derivatives `dM` where
`M=L*L'`.
# TODO make it more memory efficient
"""
function chol_gradient!{T <: Real}(g::AbstractVecOrMat{T},
  dM::AbstractVecOrMat{T}, L::AbstractMatrix{T})
  n = size(L, 1)
  At_mul_B!(g, spduplication(n),
    kron(L', speye(n)) * (dM + spcommutation(n) * dM))
end
function chol_gradient{T <: Real}(dM::AbstractVecOrMat{T}, L::AbstractMatrix{T})
  n = size(L, 1)
  if ndims(dM) == 1 # vector
    g = zeros(T, binomial(n + 1, 2))
  else # matrix
    g = zeros(T, binomial(n + 1, 2), size(dM, 2))
  end
  chol_gradient!(g, dM, L)
end

"""
    kronaxpy!(A, X, Y)

Overwrites `Y` with `A ⊗ X + Y`. Same as `Y += kron(A, X)` but more efficient.
"""
function kronaxpy!{T <: Real}(A::AbstractMatrix{T},
  X::AbstractMatrix{T}, Y::AbstractMatrix{T})

  # retrieve matrix sizes
  m, n = size(A, 1), size(A, 2)
  p, q = size(X, 1), size(X, 2)
  # loop over (i,j) blocks of Y
  irange, jrange = 1:p, 1:q
  @inbounds for j in 1:n, i in 1:m
    a = A[i, j]
    irange = ((i - 1) * p + 1):(i * p)
    jrange = ((j - 1) * q + 1):(j * q)
    Yij = sub(Y, irange, jrange)  # view of (i, j)-block
    @simd for k in eachindex(Yij)
      Yij[k] += a * X[k]
    end
  end
  Y
end

# """
# Construct an `n x n` lower triangular matrix from data `l`.
# """
# function Base.LowerTriangular(l::AbstractVecOrMat, n::Integer)
#   data = zeros(eltype(l), n, n)
#   idx = 1
#   @inbounds for j in 1:n
#     @simd for i in j:n
#       data[i, j] = l[idx]
#       idx += 1
#     end
#   end
#   LowerTriangular(data)
# end # function Base.LowerTriangular

# """
# Copy lower triangular part of `A + A' - diag(A)` into `data`.
# """
# function duplicate!(data::AbstractVecOrMat, A::AbstractMatrix)
#   n = size(A, 1)
#   @assert n == size(A, 2) "A has to be a square matrix"
#   idx = 1
#   @inbounds for j in 1:n
#     data[idx] = A[j, j]
#     idx += 1
#     @simd for i in (j + 1):n
#       data[idx] = A[i, j] + A[j, i]
#       idx += 1
#     end
#   end
# end # function Base.LowerTriangular

#end # module MultivariateCalculus

using LinearAlgebra


# Compute y = A*x quickly
#
function fast_mult(Z, x)

  # TODO: Replace this slow reference implementation
  n = length(x)
  A = Matrix{Float64}(I, n, n) + Z*Z'
  return A*x

end


# Compute kappa = cond(A) quickly
#
function fast_cond(U, S, V)

  # TODO: Replace this slow reference implementation
  n = size(U,1)
  Z = U*diagm(S)*V'
  A = Matrix{Float64}(I, n, n) + Z*Z'
  return cond(A)

end


# Compute x = A\b quickly
function fast_solve(Q, R, b)

  n = length(b)
  Z = Q*R
  A = Matrix{Float64}(I, n, n) + Z*Z'
  return A\b

end


# --- Test harness ---

n = 100
k = 5
Z = rand(Float64, n, k)
A = Matrix{Float64}(I, n, n) + Z*Z'

x = rand(Float64, n)
yref = A*x
y = fast_mult(Z, x)
relerr = norm(y-yref)/norm(yref)
print("1a: Relative err $(relerr)\n")

U, S, V = svd(Z)
kappa_ref = cond(A, 2)
kappa = fast_cond(U, S, V)
relerr = abs(kappa-kappa_ref)/kappa_ref
print("1b: Relative err $(relerr)\n")

Q, R = qr(Z)
b = rand(Float64, n)
xref = A\b
x = fast_solve(Q, R, b)
relerr = norm(x-xref)/norm(xref)
print("1c: Relative err $(relerr)\n")
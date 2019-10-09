using LinearAlgebra


# Differentiate the squared M-norm of x
#
function deriv_mnorm(x, M, dx, dM)
  # Compute derivative of squared M-norm of x
  return dnormx
end


# Differentiate the solution to (I+ZZ')x = b
#
function deriv_solve(Q, R, x, b, dZ, db)
  # Compute the derivative of the solution to the system
  return dx
end


# Differentiate the Cholesky factorization of M
#
function deriv_chol(M, R, dM)
  # Compute the derivative of the Cholesky factor
  return dR
end


# -- Finite difference check for dnormx

x = rand(Float64, 10, 1)
A = rand(Float64, 10, 10)
M = A'*A
dx = rand(Float64, 10, 1)
dM = rand(Float64, 10, 10)
dM = (dM + dM')/2;

h = 1e-6
normx_p = dot(x+h*dx, (M+h*dM)*(x+h*dx))
normx_m = dot(x-h*dx, (M-h*dM)*(x-h*dx))
dnormx_fd = (normx_p - normx_m)/(2*h)

dnormx = deriv_mnorm(x, M, dx, dM)
relerr = abs(dnormx_fd-dnormx)/abs(dnormx)
println("2a: Relerr (vs finite diff): $(relerr)")


# -- Finite difference check for solve

n = 100
k = 5
Z = rand(Float64, n, k)
dZ = rand(Float64, n, k)
Q, R = qr(Z)
b = rand(Float64, n)
db = rand(Float64, n)

Zp = Z+h*dZ
Zm = Z-h*dZ
bp = b+h*db
bm = b-h*db
x = (Matrix{Float64}(I, n, n) + Z*Z')\b
xp = (Matrix{Float64}(I, n, n) + Zp*Zp')\bp
xm = (Matrix{Float64}(I, n, n) + Zm*Zm')\bm
dx_fd = (xp - xm)/(2*h)

dx = deriv_solve(Q, R, x, b, dZ, db)
relerr = norm(dx_fd-dx)/norm(dx)
println("2b: Relerr (vs finite diff): $(relerr)")


# -- Finite difference check for Cholesky

n = 10
A = rand(Float64, n, n)
M = A'*A
R = cholesky(M).U
dM = rand(Float64, n, n)
dM = (dM + dM')/2

Rp = cholesky(M+h*dM).U
Rm = cholesky(M-h*dM).U
dR_fd = (Rp-Rm)/(2*h)

dR = deriv_chol(M, R, dM)
relerr = norm(dR_fd-dR)/norm(dR)
println("2c: Relerr (vs finite diff): $(relerr)")

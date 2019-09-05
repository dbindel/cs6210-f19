# HW 1, problem 3

using LinearAlgebra
using Printf


function hw1p3_ref(u, v, x)
  n = length(u)
  A = Matrix{Float64}(I, n, n) + u * v'

  y = A*x
  z = A'*x
  d = diag(A)
  df = diag(A[n:-1:1,:])
  t = tr(A)
  c = det(A)

  return y, z, d, df, t, c
end


function hw1p3(u, v, x)
  n = length(u)

  # Problem: Replace these lines with O(n) version
  A = Matrix{Float64}(I, n, n) + u * v'
  y = A*x
  z = A'*x
  d = diag(A)
  df = diag(A[n:-1:1,:])
  t = tr(A)
  c = det(A)

  return y, z, d, df, t, c
end


# === Test script ===

# First test: Check agreement between reference computation and O(n) version

n = 500
u = rand(n)
v = rand(n)
x = rand(n)

y1, z1, d1, df1, t1, c1 = hw1p3_ref(u, v, x)
y2, z2, d2, df2, t2, c2 = hw1p3(u, v, x)

@printf("""Rel differences:
 y:  %e
 z:  %e
 d:  %e
 df: %e
 t:  %e
 c:  %e
""",
        norm(y1-y2)/norm(y1), norm(z1-z2)/norm(z1),
        norm(d1-d2)/norm(d1), norm(df1-df2)/norm(df1),
        abs(t1-t2)/abs(t1), abs(c1-c2)/abs(c1))

# Second test: Check scalability (probably don't uncomment until fixed!)

n = 500000
u = rand(n)
v = rand(n)
x = rand(n)

#y2, z2, d2, df2, t2, c2 = hw1p3(u, v, x)

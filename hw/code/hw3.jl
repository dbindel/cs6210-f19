# CS 6210, HW 3

using LinearAlgebra
using Printf


# --- Your implementations (replace reference implementations) ---


# Compute the two-norm condition number of the Gauss transform
#   G = I - tau * ek'
# where tau(1:k) = 0 and ek denotes the kth column of the identity.
# This implementation should run in O(n) time and be numerically
# stable in the face of large or small tau.
#
function hw3p1(k, tau)
  n = length(tau)
  G = Matrix{Float64}(I, n, n)
  G[k+1:end,k] = tau[k+1:end]
  kappa = cond(G)
  return kappa
end


# Compute the determinant of A and the solution to Ax = y for
#   A = [diag(d) b; c' f]
# This version of the computation should run in O(n) time.
#
function hw3p2(d, b, c, f, y)
  A = [diagm(d) b ; c' f]
  detA = det(A)
  x = A\y
  return detA, x
end


# --- Reference implementations ---


# Compute the two-norm condition number of the Gauss transform
#   G = I - tau * ek'
# where tau(1:k) = 0 and ek denotes the kth column of the identity.
#
function hw3p1_ref(k, tau)
  n = length(tau)
  G = Matrix{Float64}(I, n, n)
  G[k+1:end,k] = tau[k+1:end]
  kappa = cond(G)
  return kappa
end


# Compute the determinant of A and the solution to Ax = y for
#   A = [diag(d) b; c' f]
#
function hw3p2_ref(d, b, c, f, y)
  A = [diagm(d) b ; c' f]
  detA = det(A)
  x = A\y
  return detA, x
end


# --- Basic testers ---


function hw3p1_test()
  tau = zeros(10,1)
  kappa_ref = 1
  kappa     = hw3p1(3, tau)
  @printf("Rel err when tau = 0: %e\n", abs(kappa_ref-kappa)/abs(kappa));

  tau[4:10] = randn(7);
  kappa_ref = hw3p1_ref(3, tau);
  kappa     = hw3p1(3, tau);
  @printf("Rel err when tau = rand: %e\n", abs(kappa_ref-kappa)/abs(kappa));

  # TODO: It may make sense to add a test here for tau large (~10^16)
  #       This may require some care, as the reference calculation is not
  #       guaranteed to obtain high relative accuracy in this case!

end


function hw3p2_test()
  n = 10;
  d = rand(n);
  b = rand(n);
  c = rand(n);
  f = rand();
  y = rand(n+1);
  detA_ref, x_ref = hw3p2_ref(d, b, c, f, y);
  detA, x = hw3p2(d, b, c, f, y);

  @printf("relerr detA: %e\n", abs(detA_ref-detA)/abs(detA));
  @printf("relerr x:    %e\n", norm(x-x_ref)/norm(x));
end


hw3p1_test()
hw3p2_test()

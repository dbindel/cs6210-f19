using LinearAlgebra
using Printf

function hw7newton(A, b, c, sigma)

  n = length(b)
  In = Matrix{Float64}(I, n, n)
  en = zeros(n+1)
  en[end] = 1

  gs = zeros(4,1);
  for k = 1:4

    M = [A-sigma*In b ;
         c' 0];
    F = factorize(M)

    # Save g(sigma) into gs[k] and update sigma to new
    # eigenvalue estimate.

  end

  return sigma, gs
end


A = rand(10,10)
b = rand(10)
c = rand(10)
sigma = 5

sigma, gs = hw7newton(A, b, c, sigma)
err = minimum(abs.(eigvals(A).-sigma))
@printf("Final error: %e\n", err)
@printf("Values for g:\n")
for gk in gs
  @printf("  %e\n", gk)
end
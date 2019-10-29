using LinearAlgebra
using Printf

# Compute the Hessenberg decomposition of H+uv' in O(n^2)
#
function hw7hess(H0, u, v)

  # Replace this with something more efficient!
  F = hessenberg(H0 + u*v')
  return F.Q, F.H

end


H0 = triu(rand(10,10),-1);
u = rand(10);
v = rand(10);
A = H0+u*v';

Q, H = hw7hess(H0, u, v)
err1 = norm(A-Q*H*Q')/norm(A)
err2 = norm(tril(H,-2))
err3 = norm(Q'*Q-Matrix{Float64}(I,10,10))

@printf("Check recovery of A: %e\n", err1)
@printf("Check Hessenberg H:  %e\n", err2)
@printf("Check orthogonal Q:  %e\n", err3)

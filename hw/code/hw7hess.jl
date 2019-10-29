using LinearAlgebra
using Printf

# Compute x = (H-sigma*I)\b in O(n^2) time
#
function hw7hess(H, sigma, b)

  # Replace this with something more efficient!
  In = Matrix{Float64}(I, length(b), length(b))
  return (H-sigma*In)\b

end


H = triu(rand(10,10),-1);
b = rand(10);
sigma = rand();

In = Matrix{Float64}(I, 10, 10)
xref = (H-sigma*In)\b
x = hw7hess(H, sigma, b)

@printf("Relerr: %e\n", norm(x-xref)/norm(xref))

using LinearAlgebra


# Julia does not seem to have a direct (built-in) analogue of polyeig,
# so we will build our own based on a simple linearization.  We will
# keep only the real eigenvalues (the complex ones have no meaning in
# this context).
#
function multiplier_qep(A, b)
  n = length(b)
  C = [zeros(n, n)  Matrix{Float64}(I, n, n) ;
       b*b'-A*A     2*A]
  mu = eigvals(C)
  return sort(real(mu[imag(mu) .== 0]))
end


# Return a global minimizer for the quadratically constrained
# quadratic program
#   minimize x'Ax/2 - x'b s.t. x'x = 1
# Your code should run in an overall time of O(n^3).
#
function qcqp_solve(A, b)
   # TODO
end


# Sanity checks -- I recommend checking that
#  a.  The multipliers correspond to stationary points
#  b.  The solution returned by qcqp_solve is a constrained minimizer
#  c.  For some reference problem, qcqp_solve returns the global min

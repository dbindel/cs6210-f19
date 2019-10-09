% Compute the norm of the tridiagonal inverse
%
function [normTinv] = p3_norminv(alpha, beta)

  T = diag(alpha) + diag(beta, 1) + diag(beta, -1);
  normTinv = norm(inv(T), 1);

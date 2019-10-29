% Solve the shifted Hessenberg system (H-s*I)*x = b in O(n^2)
%
function [x] = hw7hess(H, sigma, b)

  % Replace this with something more efficient!
  n = length(b);
  x = (H-sigma*eye(n))\b;

% Compute the Hessenberg decomposition of R+uv' in O(n^2)
%
function [Q, H] = hw7hess(H0, u, v)

  % Replace this with something more efficient!
  [Q, H] = hess(H0 + u*v');

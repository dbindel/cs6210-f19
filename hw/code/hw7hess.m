% Compute the Hessenberg decomposition of R+uv' in O(n^2)
%
function [Q, H] = hw7hess(R, u, v)

  % Replace this with something more efficient!
  [Q, H] = hess(R + u*v');

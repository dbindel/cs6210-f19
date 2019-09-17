% kappa = hw3p1(k, tau)
%
% Compute the two-norm condition number of the Gauss transform
%   G = I - tau*ek'
% where tau(1:k) = 0 and ek denotes the kth column of the identity.
% This implementation should run in O(n) time and be numerically
% stable in the face of large or small tau.
%
function kappa = hw3p1(k, tau)

  % Replace this code with something better!
  G = eye(length(tau));
  G(k+1:end,k) = -tau(k+1:end);
  kappa = cond(G);

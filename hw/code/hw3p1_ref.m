% kappa = hw3p1_ref(k, tau)
%
% Compute the two-norm condition number of the Gauss transform
%   G = I - tau*ek'
% where tau(1:k) = 0 and ek denotes the kth column of the identity.
%
function kappa = hw3p1_ref(k, tau)

  G = eye(length(tau));
  G(k+1:end,k) = -tau(k+1:end);
  kappa = cond(G);

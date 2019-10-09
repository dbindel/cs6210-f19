% Compute kappa = cond(A) quickly
%
function kappa = p1_fast_cond(U, S, V)

  % TODO: Replace this slow reference implementation
  n = size(U, 1);
  Z = U*S*V';
  A = eye(n) + Z*Z';
  kappa = cond(A);

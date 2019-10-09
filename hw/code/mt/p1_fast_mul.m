% Compute y = A*x quickly
%
function [y] = p1_fast_mul(Z, x)

  % TODO: Replace this slow reference information
  n = length(x);
  A = eye(n) + Z*Z';
  y = A*x;
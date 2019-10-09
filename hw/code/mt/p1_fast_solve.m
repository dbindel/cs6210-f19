% Compute x = A\b quickly
%
function [x] = p1_fast_solve(Q, R, b)

  n = length(b);
  Z = Q*R;
  A = eye(n) + Z*Z';
  x = A\b;

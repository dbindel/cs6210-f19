% [detA, x] = hw3p2(d, b, c, f, y)
%
% Compute the determinant of A and the solution to Ax = y for
%    A = [diag(d), b; c', f]
% This version of the computation should run in O(n) time.
%
function [detA, x] = hw3p2(d, b, c, f, y)

  % Replace this with a fast code!
  A = [diag(d), b; c', f];
  detA = det(A);
  x = A\y;

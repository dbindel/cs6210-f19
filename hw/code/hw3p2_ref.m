% [detA, x] = hw3p2_ref(d, b, c, f, y)
%
% Compute the determinant of A and the solution to Ax = y for
%    A = [diag(d), b; c', f]
%
function [detA, x] = hw3p2_ref(d, b, c, f, y)

  A = [diag(d), b; c', f];
  detA = det(A);
  x = A\y;

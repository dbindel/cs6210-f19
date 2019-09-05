function [y, z, d, df, t, c] = hw1p3(u, v, x)
  n = length(u);

  % Problem: replace these lines with O(n) version
  A = eye(n) + u*v';
  y  = A*x;
  z  = A'*x;
  d  = diag(A);
  df = diag(flipud(A));
  t  = trace(A);
  c  = det(A);

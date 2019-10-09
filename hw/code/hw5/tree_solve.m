% [x] = tree_solve(p, d, w, b)
%
% Solve Ax = b where A is a symmetric positive definite matrix
% with diagonal d and off-diagonal entries A(i,p(i)) = A(p(i),i) = w(i)
% (except where p(i) = 0).  Your solution should run in O(n) time, and
% should *not* use the sparse solvers in MATLAB or Julia directly --
% write your own tree factorization!  You may assume that p is in ascending
% order (i.e. the nodes are ordered from leaves to root).
%
function [x] = tree_solve(p, d, w, b)

  A = diag(d);
  for i = 1:length(p)
    if p(i) > 0
      A(i, p(i)) = w(i);
      A(p(i), i) = w(i);
    end
  end
  x = A\b;

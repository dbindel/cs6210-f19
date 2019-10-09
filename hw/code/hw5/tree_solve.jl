# [x] = tree_solve(p, d, w, b)
#
# Solve Ax = b where A is a symmetric positive definite matrix
# with diagonal d and off-diagonal entries A(i,p(i)) = A(p(i),i) = w(i)
# (except where p(i) = 0).  Your solution should run in O(n) time, and
# should *not* use the sparse solvers in MATLAB or Julia directly --
# write your own tree factorization!  You may assume that p is in ascending
# order (i.e. the nodes are ordered from leaves to root).
#
function tree_solve(p, d, w, b)
  A = zeros(Float64, length(p), length(p))
  for i = 1:length(p)
    A[i, i] = d[i]
    if p[i] > 0
      A[i, p[i]] = w[i]
      A[p[i], i] = w[i]
    end
  end
  return A\b
eend

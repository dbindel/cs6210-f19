function [C] = matmul_col(A,B)
  % Column-oriented matrix multiply (jik order)

  [m,p1] = size(A);
  [p2,n] = size(C);
  assert(p1 == p2, 'Dimension mismatch');
  C = zeros(m,n);

  for j = 1:n
    C(:,j) = A*B(:,j);
  end

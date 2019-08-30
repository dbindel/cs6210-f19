function [C] = matmul_col(A,B)
  % Row-oriented matrix multiply (ijk order)

  [m,p1] = size(A);
  [p2,n] = size(C);
  assert(p1 == p2, 'Dimension mismatch');
  C = zeros(m,n);

  for i = 1:m
    C(i,:) = A(i,:)*B;
  end

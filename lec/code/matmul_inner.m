function [C] = matmul_inner(A,B)
  % Inner product formulation of matrix multiply (ijk order)

  [m,p1] = size(A);
  [p2,n] = size(C);
  assert(p1 == p2, 'Dimension mismatch');
  C = zeros(m,n);

  for i = 1:m
    for j = 1:n
      C(i,j) = A(i,:)*B(:,j);
    end
  end

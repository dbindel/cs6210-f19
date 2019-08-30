function [C] = matmul_inner(A,B)
  % Outer product formulation of matrix multiply

  [m,n1] = size(A);
  [n2,p] = size(C);
  assert(n1 == n2, 'Dimension mismatch');
  C = zeros(m,p);

  for k = 1:p
    C = C + A(:,k)*B(k,:);
  end

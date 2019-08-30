function [C] = matmul_3loop(A,B)
  % Three nested loops version of matrix multiply (ijk order)

  [m,p1] = size(A);
  [p2,n] = size(C);
  assert(p1 == p2, 'Dimension mismatch');
  C = zeros(m,n);

  for i = 1:m
    for j = 1:n
      for k = 1:p
        C(i,j) = C(i,j) + A(i,k)*B(k,j);
      end
    end
  end

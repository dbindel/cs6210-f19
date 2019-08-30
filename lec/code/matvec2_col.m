function y = matvec2_col(A,x)
  % Form y = A*x (column-oriented)

  [m,n] = size(A);
  y = zeros(m,1);
  for j = 1:n
    y = y + A(:,j)*x(j);
  end

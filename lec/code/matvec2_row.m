function y = matvec2_row(A,x)
  % Form y = A*x (row-oriented)

  [m,n] = size(A);
  y = zeros(m,1);
  for i = 1:m
    y(i) = A(i,:)*x;
  end

function [result] = mydot(x,y)
  n = length(x);
  result = 0;
  for i = 1:n
    result = result + x(i)*y(i); % two flops/iteration
  end

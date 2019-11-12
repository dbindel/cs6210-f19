% [U] = sweep_gs(U, F)
%
% Run one Gauss-Seidel sweep for 2D Poisson
%
function [U] = sweep_gs(U, F);

  n  = size(U,1)-2;
  h2 = 1/(n+1)^2;

  for j = 2:n+1
    for i = 2:n+1
      U(i,j) = (U(i-1,j) + U(i+1,j) + ...
                U(i,j-1) + U(i,j+1) + h2*F(i,j))/4;
    end
  end

end

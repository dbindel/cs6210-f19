% [Unew] = sweep_jacobi(U, F)
%
% Run one Jacobi sweep for 2D Poisson
%
function [Un] = sweep_jacobi(U, F);

  n  = size(U,1)-2;
  h2 = 1/(n+1)^2;
  Un = U;

  for j = 2:n+1
    for i = 2:n+1
      Un(i,j) = (U(i-1,j) + U(i+1,j) + ...
                 U(i,j-1) + U(i,j+1) + h2*F(i,j))/4;
    end
  end

end

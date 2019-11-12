% [U] = sweep_gs_rb(U, F)
%
% Run one Gauss-Seidel sweep for 2D Poisson (red-black ordering)
%
function [U] = sweep_gs_rb(U, F);

  n  = size(U,1)-2;
  h2 = 1/(n+1)^2;

  for c = 0:1
    for j = 2:n+1
      for i = 2:n+1
        if mod(i+j,2) == c
          U(i,j) = (U(i-1,j) + U(i+1,j) + ...
                    U(i,j-1) + U(i,j+1) + h2*F(i,j))/4;
        end
      end
    end
  end

end

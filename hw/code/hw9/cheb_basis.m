% [V] = cheb_basis(n, cs)
% 
% Form a cs-column matrix of Chebyshev polynomials sampled on a uniform
% grid of n points.
%
function [V] = cheb_basis(n, cs)
  V = zeros(n, cs);
  xx = linspace(-1,1,n)';
  V(:,1) = 1;
  V(:,2) = xx;
  for j = 3:cs
    V(:,j) = 2*(xx.*V(:,j-1)) - V(:,j-2);
  end
end

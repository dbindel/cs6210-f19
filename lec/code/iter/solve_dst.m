% [X] = solve_dst(F)
%
% Use discrete sine transforms to solve the 2D Poisson problem
%   T*X + X*T = F
%
function [X] = solve_dst(F)
  n = size(F,1);
  lambda = zeros(n,1);
  lambda(:) = 2*(1-cos(pi/(n+1)*(1:n)));
  e = ones(n,1);
  L = lambda*e' + e*lambda';
  X = idst2d(dst2d(F)./L);
end

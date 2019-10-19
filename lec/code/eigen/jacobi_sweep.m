% [A] = jacobi_sweep(A, nsweeps)
%
% Apply nsweeps Jacobi iteration sweeps (default is 1).
%
function [A] = jacobi_sweep(A, nsweeps)

  n = length(A);
  if nargin < 2, nsweeps = 1; end

  for sweep = 1:nsweeps
    for k = 2:n
      for l = 1:k-1
        J = jacobi_rot(A(k,k), A(k,l), A(l,l));
        A([k l], :) = J'*A([k l], :);
        A(:, [k l]) = A(:, [k l])*J;
      end
    end
  end

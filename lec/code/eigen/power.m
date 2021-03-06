% [v,lambda] = power(A, v, maxiter, rtol)
%
% Run power iteration to compute the dominant eigenvalue of A and
% an associated eigenvector.  This will fail in general if there are
% multiple dominant eigenvalues (e.g. from a complex conjugate pair).
%
% Inputs:
%   A: Matrix to be analyzed
%   v: Start vector (default random)
%   maxiter: Maximum number of iterations allowed (default 1000)
%   rtol: Rel residual tolerance for convergence (default 1e-3)
%
function [v,lambda] = power(A, v, maxiter, rtol)

  % Fill in default parameters
  if nargin < 2, v = [];         end
  if nargin < 3, maxiter = 1000; end
  if nargin < 4, rtol = 1e-3;    end

  % Start with a random vector by default
  if isempty(v)
    v = randn(length(A),1);
  end
  v = v / norm(v);

  % Get estimate of 2-norm of A for convergence test
  normAest = sqrt(norm(A,1) * norm(A,inf));

  % Run the iteration
  Av = A*v;
  for k = 1:maxiter

    % Take a power method step and compute Rayleigh quotient
    v = Av/norm(Av);
    Av = A*v;
    lambda = v'*Av;

    % Compute the residual and check vs tolerance
    r = Av-v*lambda;
    normr = norm(r);
    if normr < rtol*normAest,
      return;
    end

  end

  % If we get here, give a warning
  warning(...
    sprintf('Power did not converge (rel resid %e after %d steps)', ...
            normr/normAest, maxiter));

end

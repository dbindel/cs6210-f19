% [v,lambda] = rqi(A, sigma, v, maxiter, rtol)
%
% Run Rayleigh quotient iteration to compute an eigenpair of A.
%
% Inputs:
%   A: Matrix to be analyzed
%   sigma: Initial shift
%   v: Start vector (default random)
%   maxiter: Maximum number of iterations allowed (default 1000)
%   rtol: Rel residual tolerance for convergence (default 1e-8)
%
function [v,lambda] = rqi(A, lambda, v, maxiter, rtol)

  % Fill in default parameters
  if nargin < 2, lambda = 0;     end
  if nargin < 3, v = [];         end
  if nargin < 4, maxiter = 1000; end
  if nargin < 5, rtol = 1e-8;    end

  % Start with a random vector by default
  if isempty(v)
    v = randn(length(A),1);
  end
  v = v / norm(v);

  % Get estimate of 2-norm of A for convergence test
  normAest = sqrt(norm(A,1) * norm(A,inf));
  I = eye(length(v));

  % Run the iteration
  for k = 1:maxiter

    % Take a power method step and compute Rayleigh quotient
    v = (A-lambda*I)\v;
    v = v/norm(v);
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
    sprintf('RQI did not converge (rel resid %e after %d steps)', ...
            normr/normAest, maxiter));

end

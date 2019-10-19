% [v,lambda] = subspace(A, k, V, maxiter, rtol)
%
% Orthogonal iteration to compute the dominant invariant subspace of A.
%
% Inputs:
%   A: Matrix to be analyzed
%   k: Dimension of subspace
%   v: Start vector (default random)
%   maxiter: Maximum number of iterations allowed (default 1000)
%   rtol: Rel residual tolerance for convergence (default 1e-4)
%
function [V,L] = subspace(A, k, V, maxiter, rtol)

  % Fill in default parameters
  if nargin < 2, k = 1;          end
  if nargin < 3, V = [];         end
  if nargin < 4, maxiter = 1000; end
  if nargin < 5, rtol = 1e-4;    end

  % Start with a random vector by default
  if isempty(V)
    V = randn(length(A),k);
  end
  V = V / norm(V);

  % Get estimate of 2-norm of A for convergence test
  normAest = sqrt(norm(A,1) * norm(A,inf));

  % Run the iteration
  AV = A*V;
  for k = 1:maxiter

    % Take a power method step and compute Rayleigh quotient
    [V,R] = qr(AV,0);
    AV = A*V;
    L = V'*AV;

    % Compute the residual and check vs tolerance
    R = AV-V*L;
    normR = norm(R, 'fro');
    if normR < rtol*normAest,
      return;
    end

  end

  % If we get here, give a warning
  warning(...
    sprintf('Power did not converge (rel resid %e after %d steps)', ...
            normR/normAest, maxiter));

end

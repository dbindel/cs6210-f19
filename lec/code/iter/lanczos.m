% [Q,alpha,beta] = lanczos(A,b)
%
% Compute an Lanczos decomposition
%
%   A*Q(:,1:end-1) = Q*T
%
% where T is a k+1-by-k tridiagonal matrix with diagonal
% entries alpha and super/subdiagonals beta, and Q has
% orthonormal columns.
%
function [Q,H] = lanczos(A,b,k)

  n = length(A);
  Q = zeros(n,k+1);   % Orthonormal basis
  alpha = zeros(k,1);
  beta  = zeros(k,1);

  Q(:,1) = b/norm(b);
  for j = 1:k
    Q(:,j+1) = A*Q(:,j);
    alpha(j) = Q(:,j)'*Q(:,j+1);
    Q(:,j+1) = Q(:,j+1)-alpha(j)*Q(:,j);
    if j > 1
      Q(:,j+1) = Q(:,j+1)-beta(j-1)*Q(:,j-1);
    end
    beta(j) = norm(Q(:,j+1));
    Q(:,j+1) = Q(:,j+1)/beta(j);
  end

end

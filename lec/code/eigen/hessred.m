% [H,Q] = hessred(A)
%
% Compute the Hessenberg decomposition H = Q'*A*Q using
% Householder transformations.
%
function [H,Q] = hessred(A)

  n = length(A);
  Q = eye(n);      % Orthogonal transform so far
  H = A;           % Transformed matrix so far

  for j = 1:n-2

    % -- Find W = I-2vv' to put zeros below H(j+1,j)
    u    = H(j+1:end,j);
    u(1) = u(1) + sign(u(1))*norm(u);
    v    = u/norm(u);

    % -- H := WHW', Q := QW
    H(j+1:end,:) = H(j+1:end,:)-2*v*(v'*H(j+1:end,:));
    H(:,j+1:end) = H(:,j+1:end)-(H(:,j+1:end)*(2*v))*v';
    Q(:,j+1:end) = Q(:,j+1:end)-(Q(:,j+1:end)*(2*v))*v';

  end

end

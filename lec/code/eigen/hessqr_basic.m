% [H] = hessqr_basic(H)
%
% Compute one basic (unshifted) implicit Hessenberg QR step via
% Householder transformations.
%
function H = hessqr_basic(H)

  n = length(H);
  V = zeros(2,n-1);

  % Compute the QR factorization
  for j = 1:n-1

    % -- Find W_j = I-2vv' to put zero into H(j+1,j)
    u      = H(j:j+1,j);
    u(1)   = u(1) + sign(u(1))*norm(u);
    v      = u/norm(u);
    V(:,j) = v;

    % -- H := W_j H
    H(j:j+1,:) = H(j:j+1,:)-2*v*(v'*H(j:j+1,:));

  end

  % Compute RQ
  for j = 1:n-1

    % -- H := WHW', Q := QW
    v = V(:,j);
    H(:,j:j+1) = H(:,j:j+1)-(H(:,j:j+1)*(2*v))*v';

  end

end

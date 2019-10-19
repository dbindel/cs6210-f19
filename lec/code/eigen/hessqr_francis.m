% [H] = hessqr_francis(H)
%
% Compute a (double) implicit Hessenberg QR step with Francis shift.
% Compare to hessqr_basic.
%
function [H] = hessqr_francis(H)
  % Implicit QR step using a Francis double shift
  % (there should really be some re-scalings for floating point)

  % Compute double-shift poly and initial column of H^2 + b*H + c*I
  [b,c]   = francis_poly(H);
  C1      = H(1:3,1:2)*H(1:2,1);
  C1(1:2) = C1(1:2) + b*H(1:2,1);
  C1(1)   = C1(1)   + c;

  % Apply a similarity associated with the first step of QR on C
  v        = house(C1);
  H(1:3,:) = H(1:3,:)-2*v*(v'*H(1:3,:));
  H(:,1:3) = H(:,1:3)-(H(:,1:3)*(2*v))*v';

  % Do "bulge chasing" to return to Hessenberg form
  n = length(H);
  for j = 1:n-2
    k = min(j+3,n);

    % -- Find W = I-2vv' to put zeros below H(j+1,j), H := WHW'
    v          = house(H(j+1:k,j));
    H(j+1:k,:) = H(j+1:k,:)-2*v*(v'*H(j+1:k,:));
    H(:,j+1:k) = H(:,j+1:k)-(H(:,j+1:k)*(2*v))*v';
    H(k,j)     = 0;

  end

end

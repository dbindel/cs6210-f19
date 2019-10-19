% [H] = hessqr(H)
%
% Toy implementation of Hessenberg QR iteration with Francis double
% shift strategy and deflation.
%
function [H] = hessqr(H)

  n = length(H);
  tol = norm(H,'fro') * 1e-8;
  k = 0;
  while n > 2
    if abs(H(n,n-1)) < tol
      fprintf('At step %d: Deflated 1-by-1 block\n', k);
      H(n,n-1) = 0;
      n = n-1;
    elseif abs(H(n-1,n-2)) < tol
      fprintf('At step %d: Deflated 2-by-2 block\n', k);
      H(n-1,n-2) = 0;
      n = n-2;
    else
      H(1:n,1:n) = hessqr_francis(H(1:n,1:n));
      k = k+1;
    end
  end

end

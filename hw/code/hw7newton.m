function [sigma, gs] = hw7newton(A, b, c, sigma)

  n = length(b);
  I = eye(n);
  en = zeros(n+1,1);
  en(end) = 1;

  gs = zeros(4,1);
  for k = 1:4

    M = [A-sigma*I, b; c', 0];
    [L, U] = lu(M);

    % TODO:
    %   Compute g(sigma) and save to gs(k), then 
    %   compute a Newton step to update the
    %   eigenvalue estimate sigma.

  end

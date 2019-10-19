% [J] = jacobi(alpha, beta, gamma)
%
% Compute the Jacobi rotation
%
%   [ c -s ] [ alpha beta  ] [  c s ] = [ alpha_new 0        ]
%   [ s  c ] [ beta  gamma ] [ -s c ]   [ 0         beta_new ]
%
function [J] = jacobi_rot(alpha, beta, gamma)

  if beta == 0
    J = [1, 0; 0, 1];
  else
    b = (gamma-alpha)/2/beta;
    t = sign(b)/(abs(b) + sqrt(b^2)+1);
    c = 1/sqrt(1+t^2);
    s = c*t;
    J = [c, s; -s, c];
  end


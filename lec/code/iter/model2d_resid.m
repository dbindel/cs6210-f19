% [R] = model2d_resid(U, F)
%
% Compute residual r = T*u-h^2*f
%
function [R] = model2d_resid(U, F);

  n  = size(U,1);
  h2 = 1/(n+1)^2;
  I  = 2:n+1;

  R = 4*U(I,I) - U(I,I-1) - U(I,I+1) ...
               - U(I-1,I) - U(I+1,I) - h2*F(I,I);

end

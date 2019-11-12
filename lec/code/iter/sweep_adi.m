% [U] = sweep_adi(U, F)
%
% Run one ADI sweep for 2D Poisson (single shift)
%
function [U] = sweep_adi(U, F);
  n  = size(U,1)-2;
  I  = 2:n+1;
  h2 = 1/(n+1)^2;
  dt = 1/(n+1);
  Ts = spdiags(ones(n,1)*[-1, 2+dt, -1], [-1, 0, 1], n, n);

  % Iterate on Ts*U + U*Ts = h^2*F + 2*dt*U where Ts = T+dt*I:
  %   Ts*U = h^2*F + 2*dt*U - U*Ts
  %   U*Ts = h^2*F + 2*dt*U - Ts*U
  U(I,I) = Ts\( h2*F(I,I) - U(I,I)*Ts + 2*dt*U(I,I) );
  U(I,I) = ( h2*F(I,I) - Ts*U(I,I) + 2*dt*U(I,I) )/Ts;

end

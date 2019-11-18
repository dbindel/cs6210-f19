% plot_err(x, u, uh)
%
% Plots the error in uh as an approximation to u (both on the x grid);
% also plots the DST of the error in order to understand freq-by-freq.
%
function plot_err(x, u, uh)
  subplot(2,1,1);
  plot(x, u, x, uh);
  subplot(2,1,2);
  semilogy(abs(dst(u(2:end-1)-uh(2:end-1))));
end

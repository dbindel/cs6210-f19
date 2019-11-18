% [uh] = solve_coarse(T, uh, hf, V)
%
% Compute a correction to the guess uh to the solution of T*u = hf
% by approximating the correction equation T*du = (hf-T*u)
% via Bubnov-Galerkin with the projection space V.
%
function [uh] = solve_coarse(T, uh, hf, V)
  r  = hf-T*uh;
  % TODO: Updatte uh with a coarse correction
end


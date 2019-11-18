% [uh] = solve_bgs(T, uh, hf, bs)

% Apply block Gauss-Seidel to update a guess uh for the solution to the
% linear systetm T*u = hf.  bs is the size of the Gauss-Seidel blocks.
%
function [uh] = solve_bgs(T, uh, hf, bs)
  n = length(uh)-2;
  for j = 2:bs:n+1
    I = j:min(j+bs-1,n+1);
    % TODO: Relax the entire block of nodes indicated by I
  end
end

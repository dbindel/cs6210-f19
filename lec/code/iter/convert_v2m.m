% [U] = convert_v2m(u)
%
% Convert from a vec'd representation of a field on a 2D mesh
% to a 2D array representation (with zero padding at boundary)
%
function [U] = convert_v2m(u)
  n = sqrt(length(u));
  U = zeros(n+2,n+2);
  U(2:n+1,2:n+1) = reshape(u,n,n);
end

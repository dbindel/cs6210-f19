% [u] = convert_m2v(U)
%
% Convert from a 2D mesh representation of a field to a vec'd
% representation
%
function [u] = convert_v2m(U)
  n = size(U,1)-2;
  u = reshape(U(2:end-1,2:end-1), n*n, 1);
end

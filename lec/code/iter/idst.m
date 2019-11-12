% [TX] = idst(X)
%
% The inverse DST is X = 2/(n+1) * dst(TX).
%
function [X] = idst(TX)
  m = size(TX,1);
  X = 2/(m+1) * dst(TX);
end

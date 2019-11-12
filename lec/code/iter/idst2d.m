% [TX] = idst2d(X)
%
% The inverse DST is X = 4/(n+1)/(m+1) * dst2d(TX).
%
function [X] = idst2d(TX)
  [m,n] = size(TX);
  X = 4/prod((m+1)*(n+1)) * dst2d(TX);
end

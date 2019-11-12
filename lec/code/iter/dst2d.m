function [TX] = dst2d(X)
  TX = dst(dst(X')');
end

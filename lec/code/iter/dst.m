% [TX] = dst(X)
%
% Use FFTs to compute T*X where T is the discrete sine transform
% matrix with entries
%
%  T(i,j) = sin(i*j*pi/(n+1));
% 
% The inverse DST is X = 2/(n+1) * dst(TX).
%
function [TX] = dst(X)
  [m,n] = size(X);
  Y = zeros(2*m+2,n);
  Y(2:m+1,:) = X;
  Y(m+3:end,:) = -flipud(X);
  TX = imag(fft(Y)/2);
  TX = TX(2:m+1,:);
end

% [b,c] = francis_poly(H)
%
% Compute b, c s.t. z^2 + b*z + c = (z-sigma)(z-conj(sigma))
% where sigma is the Francis double shift for H.
%
function [b,c] = francis_poly(H)

  % Get shifts via trailing submatrix
  HH    = H(end-1:end,end-1:end);
  trHH  = HH(1,1)+HH(2,2);
  detHH = HH(1,1)*HH(2,2)-HH(1,2)*HH(2,1);

  if trHH^2 > 4*detHH   % Real eigenvalues

    % Use the one closer to H(n,n)
    lHH(1) = (trHH + sqrt(trHH^2-4*detHH))/2;
    lHH(2) = (trHH - sqrt(trHH^2-4*detHH))/2;
    if abs(lHH(1)-H(end,end)) < abs(lHH(2)-H(end,end))
      lHH(2) = lHH(1);
    else
      lHH(1) = lHH(2);
    end

    % z^2 + bz + c = (z-sigma_1)(z-sigma_2)
    b = -lHH(1)-lHH(2);
    c =  lHH(1)*lHH(2);

  else

    % In the complex case, we want the char poly for HH
    b = -trHH;
    c = detHH;

  end

function [vmax] = mymax1(v)
  % Find the maximum element of a vector -- naive about NaN

  vmax = v(1);
  for k = 2:length(v)
    if v(k) > vmax, vmax = v(k); end
  end

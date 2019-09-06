function [vmax] = mymax2(v)
  % Find the max a vector -- or NaN if any element is NaN

  vmax = v(1);
  for k = 2:length(v)
    if isnan(v(k))
      vmax = v(k);
    elseif v(k) > vmax
      vmax = v(k);
    end
  end

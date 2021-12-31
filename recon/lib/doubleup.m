function y = doubleup(x, gamma)

x = squeeze(x);
y = [x; x];

y(1:2:end, :) = x.*(1/2+gamma(1)/2) + [x(1, :); x(1:end-1, :)].*(1/4-gamma(1)/4) + ...
                [x(2:end, :); x(end,:)].*(1/4-gamma(1)/4);
  
y(2:2:end, :) = [x(1, :); x(1:end-1, :)].*(sqrt(2)/8-gamma(2)/8) + x.*(1/2-sqrt(2)/8+gamma(2)/8) + ...
                [x(2:end, :); x(end,:)].*(1/2-sqrt(2)/8+gamma(2)/8) + ...
                [x(3:end, :); x(end, :); x(end, :)].*(sqrt(2)/8-gamma(2)/8);

end
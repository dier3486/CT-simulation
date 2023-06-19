function y = doubleup(x, Gamma)

x = squeeze(x);
y = [x; x];

y(1:2:end, :) = x.*(1/2+Gamma(1)/2) + [x(1, :); x(1:end-1, :)].*(1/4-Gamma(1)/4) + ...
                [x(2:end, :); x(end,:)].*(1/4-Gamma(1)/4);
  
% y(2:2:end, :) = [x(1, :); x(1:end-1, :)].*(sqrt(2)/8-Gamma(2)/8) + x.*(1/2-sqrt(2)/8+Gamma(2)/8) + ...
%                 [x(2:end, :); x(end,:)].*(1/2-sqrt(2)/8+Gamma(2)/8) + ...
%                 [x(3:end, :); x(end, :); x(end, :)].*(sqrt(2)/8-Gamma(2)/8);
% wrong??

y(2:2:end, :) = ([x(1, :); x(1:end-1, :)] + [x(3:end, :); x(end, :); x(end, :)]).*(1/4-sqrt(2)/8-Gamma(2)/8) + ...
                (x + [x(2:end, :); x(end,:)]).*(1/4+sqrt(2)/8+Gamma(2)/8);

end
function x = doubleup2(x, Gamma)

if nargin<2
    Gamma = [0.7000    0.8854];
end

x(2:2:end, :) = [x(1, :); x(1:2:end-3, :)].*(sqrt(2)/8-Gamma(2)/8) + x(1:2:end-1, :).*(1/2-sqrt(2)/8+Gamma(2)/8) + ...
                [x(3:2:end-1, :); x(end-1,:)].*(1/2-sqrt(2)/8+Gamma(2)/8) + ...
                [x(5:2:end-1, :); x(end-1, :); x(end-1, :)].*(sqrt(2)/8-Gamma(2)/8);
x(1:2:end, :) = x(1:2:end-1, :).*(1/2+Gamma(1)/2) + [x(1, :); x(1:2:end-3, :)].*(1/4-Gamma(1)/4) + ...
                [x(3:2:end-1, :); x(end-1,:)].*(1/4-Gamma(1)/4);

end
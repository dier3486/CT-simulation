function x = doubleupbutterfly(x, Gamma)

if nargin<2
    Gamma = [0.7000    0.8854];
end

a1 = Gamma(1);
a2 = (sqrt(2)-2+Gamma(2))/(1-Gamma(1));

% (x_0 + x_1)/2
x(2:2:end, :) = (x(1:2:end-1, :) + [x(3:2:end-1, :); x(end-1, :)])./2;

% x_1
x(1:2:end-1, :) = x(1:2:end-1, :).*a1 + ...
    ([x(1, :); x(2:2:end-2, :)] + [x(2:2:end-2, :); x(end-1, :)]).*((1-a1)/2);

% x_2
x(2:2:end, :) = x(2:2:end, :) + ...
    (x(2:2:end, :) - (x(1:2:end-1, :) + [x(3:2:end-1, :); x(end, :)])./2).*a2;

end
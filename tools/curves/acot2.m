function y = acot2(x)
% my acot
% y = acot2(x)

y = acot(x) + (1-sign(x)).*(pi/2);
end
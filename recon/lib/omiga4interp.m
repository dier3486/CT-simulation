function [t_odd, t_even, gamma] = omiga4interp(x, c)
% Sobolev space linearized omiga-4-points interpoloation
% [t_odd, t_even, gamma] = omiga4interp(x, c);
% or, [t_odd, t_even] = omiga4interp(x, c);
% after that, the omiga-4-points interplation reads
% y = interp1(y0(1:2:end), t_odd)./2 + interp1(y0(2:2:end), t_even)./2 + interp1(conv(y0, [-1 2 -1], 'same'), x).*gamma./4;

x_flr = floor(x);
s_odd = mod(x_flr, 2);
alpha = x - x_flr;
beta = 1/2-sqrt(1+alpha.*(1-alpha).*4)./2;

t1 = (1+alpha-beta)./2;
t2 = (alpha+beta)./2;

t_odd = (x_flr+s_odd)./2 + t2.*s_odd + t1.*(1-s_odd);
t_even = (x_flr-s_odd)./2 + t1.*s_odd + t2.*(1-s_odd);

if nargout>2
    gamma = c(1)./sqrt(1-alpha.*(1-alpha).*c(2));
end

end

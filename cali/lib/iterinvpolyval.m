function y = iterinvpolyval(p, x, r, n)
% return y satisfy ((((y*p(1)+1)*y*p(2)+1)*y*p(3)+1)*...)y*p(n) = x
% y = iterinvpolyval(p, x, r, n);

if nargin<3
    r = [min(x(:)./p(end)).*0.8 max(x(:)./p(end)).*1.2];
end
if nargin<4
    n = 1000;
end

yy = linspace(r(1), r(2), n);
xx = iterpolyval(p, yy);
y = interp1(xx, yy, x, 'linear', 'extrap');

end
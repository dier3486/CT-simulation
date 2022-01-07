% test of theta curve fit
% run after comptonform_test4

y0 = squeeze(A4(100,100,:));
y1 = squeeze(A4(20,180,:));
x0 = theta_x(:);

% x0 = theta(:)./theta(end);


m = 1; n = 5;
t0 = ones(1, m+n);

t1 = lsqcurvefit(@(x, xdata) pptestval(x(1:m), x(m+1:end), x0, y0, xdata), t0, x0, y1);

z1 = pptestval(t1(1:m), t1(m+1:end), x0, y0, x0);

function y = ponpval(p1, p2, x)

y = polyval(p1, x)./polyval(p2, x);

end

function y = pptestval(p1, p2, x0, y0, x)

x1 = polyval(p1, x).*(1-x).*x+x;
y = interp1(x0, y0, x1, 'linear', 'extrap').*polyval(p2, x);


end


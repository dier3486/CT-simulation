x = [randn(100,1); randn(300,1).*1.5+8];

mx = mean(x);
mdx = median(x);
Nx = size(x(:),1);

lb = min(x);
ub = max(x);
x1 = (x-lb)./(ub-lb);

k = 1.01;
% u = x1 - mean(x1);
mu = mean(x1);
mdu = median(x1);

v = fzero(@(x) zfun(x1, x, k-1), mu);

dv = zfun(x1, mu, k-1);
dv = -abs(dv*k)^(1/(k-1))*sign(dv);

function r = zfun(u, x, k)

r = mean(abs(u-x).^k.*sign(u-x));
end
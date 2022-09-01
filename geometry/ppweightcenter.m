function c = ppweightcenter(cs)
% to find the weight center of a spline
% c = ppweightcenter(cs);
% where the cs is the return of sline.m, e.g.  cs = spline(x,[0 y 0]);

dx = diff(cs.breaks);
dx = dx(:);
X = [dx.^5./5 dx.^4./4 dx.^3./3 dx.^2./2 dx];

w = sum(reshape(cs.coefs.*X(:, 2:end), [], 1));
c = cs.breaks(1:end-1)*sum(cs.coefs.*X(:, 2:end),2) + sum(reshape(cs.coefs.*X(:, 1:end-1), [], 1));
c = c/w;
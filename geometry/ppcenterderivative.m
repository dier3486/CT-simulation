function [dp, p] = ppcenterderivative(cs)
% partial derivative the weight center of a spline by the node position
% [dp, p] = ppcenterderivative(cs);
% where the cs is the return of sline.m, e.g.  cs = spline(x,[0 y 0]);
% the vectorized version of this function will be used in inverse geometry calibration

dx = diff(cs.breaks);
dx = dx(:);
X = [dx.^5./5 dx.^4./4 dx.^3./3 dx.^2./2 dx];

wi = sum(cs.coefs.*X(:, 2:end), 2);
ui = sum(cs.coefs.*X(:, 1:end-1), 2);

w = sum(wi);
c = cs.breaks(1:end-1)*wi + sum(ui);

dwdx = -diff([0; wi./dx; 0]);
dcdx = [wi; 0] - diff([0; wi.*cs.breaks(1:end-1)'./dx; 0]) - diff([0; ui./dx; 0]).*2;

dp = (dcdx.*w - dwdx.*c)./w^2;
p = c/w;

% w = sum(reshape(cs.coefs.*X(:, 2:end), [], 1));
% c = cs.breaks(1:end-1)*sum(cs.coefs.*X(:, 2:end),2) + sum(reshape(cs.coefs.*X(:, 1:end-1), [], 1));
% c = c/w;

end
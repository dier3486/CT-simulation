function [u, v] = square2circle(x, y, flag_polar)
% map points from |r|<=1 to |r|^2<=1 of 2-D.
% function [u, v] = square2circle(x, y)

if nargin<3
    flag_polar = false;
end

% rho
rho = max(abs(x), abs(y));
% theta
s1 = y>x;
s2 = y>-x;
xony = x./y;
yonx = y./x;
theta = fillmissing(yonx.*(~s1&s2), 'constant', 0);
theta = theta + fillmissing((4+yonx).*(s1&~s2), 'constant', 0);
theta = theta + fillmissing((2-xony).*(s1&s2), 'constant', 0);
theta = theta + fillmissing((6-xony).*(~s1&~s2), 'constant', 0);
theta = theta./8;

if flag_polar
    u = rho;
    v = mod(theta, 1);
else
    theta = theta.*(pi*2);
    u = rho.*cos(theta);
    v = rho.*sin(theta);
end

end
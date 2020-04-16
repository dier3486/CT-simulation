1;



function x = iterinvpolyval_t1(p, y, x0, tol)

if nargin<3
    x0 = y;
end
if nargin<4
    tol = 1e-8;
end
[m, n] = size(y);

Nmax = 100;
x = x0;
alpha = 1.0;
s = true(m, n);
% s0 = s;
for iter = 1:Nmax
    r = iterpolyval_s(p, x, s) - y(s);
    s_r = abs(r)>tol;
    s(s) = s_r;
    r = r(s_r);
    if all(~s_r(:))
        break;
    end
    dy = iterpolydiff_t1(p, x, s);
    x(s) = x(s) - r./dy.*alpha;
end
% iter
end

function dy = iterpolydiff_t1(p, x, s)

n = size(p, 2);
dy = ones(size(x));
for ii = 1:n-1
    dy(~s) = 0;
    dy(s) = dy(s).*x(s);
    dy = dy*p(:, ii).*((n-ii+1)/(n-ii)) + 1.0;
end
dy(~s) = 0;
dy = dy.*p(:, n);
dy = dy(s);

end

function dy = iterpolydiff_t0(p, x)

n = size(p, 2);
dy = ones(size(x));
for ii = 1:n-1
    dy = (dy.*x).*p(:, ii).*((n-ii+1)/(n-ii)) + 1.0;
end
dy = dy.*p(:, n);

end

function y = iterpolyval_s(p, x, s)
% return y = ((((x*p(1)+1)*x*p(2)+1)*x*p(3)+1)*...)*x*p(n)
% y = iterpolyval(p, x);

n = size(p, 2);
y = ones(size(x));
for ii = 1:n-1
    y(~s) = 0;
    y(s) = y(s).*x(s);
    y = y.*p(:, ii) + 1.0;
end
y(~s) = 0;
y(s) = y(s).*x(s);
y = y.*p(:, n);
y = y(s);
end
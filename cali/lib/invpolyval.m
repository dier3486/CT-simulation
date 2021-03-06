function x = invpolyval(p, y, x0, lrb, tol, Nmax)
% find x: y = polyval(p, x);
% x = invpolyval(p, y);

if nargin<3 || isempty(x0)
    x0 = (y-p(end))./p(end-1);
    x0 = fillmissing(x0, 'constant', 0);
end
if nargin<4
    lrb = [];
end
if nargin<5
    tol = 1e-6;
end
if nargin<6
    Nmax = 100;
end


% debug
% xc = fzero(@(x) polyval(p, x)-y, x0);

p = p(:)';
Np = length(p);
dp = (Np-1:-1:1).*p(1:end-1);

alpha = 1.0;

r0 = fillmissing((y-polyval(p, x0))./polyval(dp, x0), 'constant', 0);
x = x0;
r = r0;

for ii = 1:Nmax
    x = x + r.*alpha;
    rr = r.^2;
    if all(rr<tol^2)
        break;
    end
    r0 = r;
    r = fillmissing((y - polyval(p, x))./polyval(dp, x), 'constant', 0);
    s = sign(r)~=sign(r0);
%     if any(s)
%         1;
%     end
    r(s) = min(abs(r(s)), abs(r0(s))./2).*sign(r(s));
end

if ~isempty(lrb)
    ylr = polyval(p, lrb);
    dlr = polyval(dp, lrb);
    s1 = y<ylr(1);
    x(s1) = lrb(1) + (y(s1)-ylr(1))./dlr(1);
    s2 = y>ylr(2);
    x(s2) = lrb(2) + (y(s2)-ylr(2))./dlr(2);
end

end
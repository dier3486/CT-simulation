function y = iterpolyval(p, x)
% return y = ((((x*p(1)+1)*x*p(2)+1)*x*p(3)+1)*...)*x*p(n)
% y = iterpolyval(p, x);

n = size(p, 2);
y = ones(size(x));
for ii = 1:n-1
    y = (y.*x).*p(:, ii) + 1.0;
end
y = (y.*x).*p(:, n);

end
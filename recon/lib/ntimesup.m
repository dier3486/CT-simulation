function y = ntimesup(x, n, Gamma)

x = squeeze(x);
y = repmat(x.*0, 1, 1, n);

n = cast(n, 'like', Gamma);
a = (0:n-1)./n;
b = 1/2-sqrt(1+4.*a.*(1-a))./2;
c = Gamma(1)./sqrt(1-a.*(1-a).*Gamma(2));
t = [(2-a-b+(2-a.*3).*c)./4 (b+(1-a).*(1-c))./4];

u = reshape(x(:)*t, [size(x) n*2]);
y(:,:, 1) = u(:,:, 1) + [u(1, :, n+1); u(1:end-1, :, n+1)] + [u(2:end, :, n+1); u(end, :, n+1)];

for ii = 2:n
    y(:,:, ii) = [u(1, :, ii+n); u(1:end-1, :, ii+n)] + u(:,:, ii) + [u(2:end, :, n+2-ii); u(end, :, n+2-ii)] + ...
                 [u(3:end, :, n*2+2-ii); u(end, :, n*2+2-ii); u(end, :, n*2+2-ii)];
end

y = reshape(permute(y, [3 1 2]), size(x).*[n 1]);

end
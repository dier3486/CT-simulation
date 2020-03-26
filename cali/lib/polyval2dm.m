function Z = polyval2dm(P, X, Y)
% 2D polyvalm
% Z = polyval2dm(P, X, Y);
% Z = P(1,1) + P(2,1)*X + P(2,2)*X*Y + P(1,2)*Y + ... + P(m+1,n+1)*X^m*Y^n.

Z = zeros(size(X));

[m,n] = size(P);
[im, in] = ndgrid(1:m, 1:n);
sP = P~=0;
im = im(sP); in = in(sP);
index = im + (in-1).*m;
Ns = sum(sP(:));

for ii = 1:Ns
    Z = Z + P(index(ii)).*X.^(im(ii)-1).*Y.^(in(ii)-1);
end

end
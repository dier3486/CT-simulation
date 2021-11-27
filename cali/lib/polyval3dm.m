function R = polyval3dm(P, X, Y, Z)
% 3D polyvalm
% R = polyval3dm(P, X, Y, Z);
% R = P(1,1,1) + P(2,1,1)*X + P(1,2,1)*Y + P(1,1,2)*Z + ... + P(m+1,n+1,q+1)*X^m*Y^n*Z^q.

% R = zeros(size(X), class(X));
R = X;  R(:) = 0;

[m,n,q] = size(P);
[im, in, iq] = ndgrid(1:m, 1:n, 1:q);
sP = P~=0;
im = im(sP); in = in(sP); iq = iq(sP);
index = im + (in-1).*m + (iq-1).*(m*n);
Ns = sum(sP(:));

for ii = 1:Ns
    R = R + P(index(ii)).*X.^(im(ii)-1).*Y.^(in(ii)-1).*Z.^(iq(ii)-1);
end

end
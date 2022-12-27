function R = polyval2dm(P, X, Y)
% 2D polyvalm
% R = polyval2dm(P, X, Y);
% R = P(1,1) + P(2,1)*X + P(2,2)*X*Y + P(1,2)*Y + ... + P(m+1,n+1)*X^m*Y^n.

[m, n] = size(P);

R = zeros(size(X), 'like', X);
for ix = m:-1:1
    yy = P(ix, n);
    for iy = n-1:-1:1
        yy = yy.*Y + P(ix, iy);
    end
    R = R.*X + yy;
end

end
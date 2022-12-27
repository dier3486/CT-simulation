function R = polyval3dm(P, X, Y, Z)
% 3D polyvalm
% R = polyval3dm(P, X, Y, Z);
% R = P(1,1,1) + P(2,1,1)*X + P(1,2,1)*Y + P(1,1,2)*Z + ... + P(m+1,n+1,q+1)*X^m*Y^n*Z^q.

[m, n, q] = size(P);
R = zeros(size(X), 'like', X);
yy = zeros(size(X), 'like', X);
for ix = m:-1:1
    yy = yy.*0;
    for iy = n:-1:1
        zz = P(ix, iy, q);
        for iz = q-1:-1:1
            zz = zz.*Z + P(ix, iy, iz);
        end
        yy = yy.*Y + zz;
    end
    R = R.*X + yy;
end

end
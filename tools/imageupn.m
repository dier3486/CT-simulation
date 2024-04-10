function img2 = imageupn(img, n)
% image n-times upsampling by linear interps

[Nx, Ny, Nimg] = size(img);

% fill edges
img1 = zeros(Nx+2, Ny+2, Nimg, 'like', img);
img1(2:end-1,2:end-1,:) = img;
img1(1, 2:end-1, :) = img(1, :, :);
img1(end, 2:end-1, :) = img(end, :, :);
img1(:, 1, :) = img1(:, 2, :);
img1(:, end, :) = img1(:, end-1, :);

modd = mod(n,2);

% A0 = ((1:n)'.*2-1)*((1:n).*2-1)./n^2./4;

d = n - (1:n) + (1+modd)/2;
A = d'*d./n^2;

A = A(triu(true(n)));
A = [A' 0];

Aidx0 = zeros(n, n);
Aidx0(triu(true(n))) = 1 : n*(n+1)/2;
Aidx0 = Aidx0 + tril(Aidx0', -1);
Aidx = ones(n+1, n+1).*(n*(n+1)/2+1);
Aidx(1:n, 1:n) = Aidx0;

index_ij = zeros(1, 4);

img1_ext = reshape(img1(:)*A, Nx+2, Ny+2, Nimg, n*(n+1)/2+1);


K = circshift(1:n, floor(n/2));

img2 = zeros(Nx*n, Ny*n, Nimg);
for ii = 1:n
    for jj = 1:n
        index_ij(1) = Aidx(K(ii), K(jj));
        index_ij(2) = Aidx(n+modd-K(ii)+1, K(jj));
        index_ij(3) = Aidx(K(ii), n+modd-K(jj)+1);
        index_ij(4) = Aidx(n+modd-K(ii)+1, n+modd-K(jj)+1);
        mbase = [ii, jj] > n/2;
        img2(ii:n:end, jj:n:end, :) = ...
            img1_ext(1+mbase(1):end-2+mbase(1), 1+mbase(2):end-2+mbase(2), :, index_ij(1)) + ...
            img1_ext(2+mbase(1):end-1+mbase(1), 1+mbase(2):end-2+mbase(2), :, index_ij(2)) + ...
            img1_ext(1+mbase(1):end-2+mbase(1), 2+mbase(2):end-1+mbase(2), :, index_ij(3)) + ...
            img1_ext(2+mbase(1):end-1+mbase(1), 2+mbase(2):end-1+mbase(2), :, index_ij(4));
    end
end


end
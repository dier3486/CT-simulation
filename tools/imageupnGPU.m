function img2 = imageupnGPU(img, n)
% image n-times upsampling by linear interps

[Ny, Nx, Nimg] = size(img);

% imgClass = classGPU(img);

% fill edges
img1 = zeros(Ny+2, Nx+2, Nimg, 'like', img);
img1(2:end-1,2:end-1,:) = img;
img1(1, 2:end-1, :) = img(1, :, :);
img1(end, 2:end-1, :) = img(end, :, :);
img1(:, 1, :) = img1(:, 2, :);
img1(:, end, :) = img1(:, end-1, :);

Xgrid = ((1:Nx*n) + (3/2*n-1/2))./n;
Ygrid = ((1:Ny*n) + (3/2*n-1/2))./n;

% Xgrid = gpuArray(cast(Xgrid, imgClass));
% Ygrid = gpuArray(cast(Ygrid, imgClass));
Xgrid = cast(Xgrid, 'like', real(img(1)));
Ygrid = cast(Ygrid, 'like', real(img(1)));

[X, Y] = meshgrid(Xgrid, Ygrid);

img2 = zeros(Ny*n, Nx*n, Nimg, 'like', img);
for ii = 1:Nimg
    img2(:,:, ii) = interp2(img1(:,:, ii), X, Y);
end


end
function img2 = imagedown(img, n)
% image n-times down-sampling by mean of the n*n blocks

if nargin < 2
    n = 2;
end

[Nx, Ny, ~] = size(img);

Nx2 = ceil(Nx/n);
Ny2 = ceil(Ny/n);

if mod(Nx, n) ~= 0
    img = cat(1, img, repmat(img(end, :, :), n-mod(Nx, n), 1, 1));
end
if mod(Ny, n) ~= 0
    img = cat(2, img, repmat(img(:, end, :), 1, n-mod(Ny, n), 1));
end

img = permute(reshape(img, n, Nx2, n, Ny2, []), [1 3 2 4 5]);

img2 = squeeze(mean(reshape(img, n^2, Nx2, Ny2, []), 1));

end
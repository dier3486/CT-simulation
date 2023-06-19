function img2 = imageup2(img)

[Nx, Ny, Nimg] = size(img);

% fill edges
img1 = zeros(Nx+2, Ny+2, Nimg);
img1(2:end-1,2:end-1,:) = img;
img1(1, 2:end-1, :) = img(1, :, :);
img1(end, 2:end-1, :) = img(end, :, :);
img1(:, 1, :) = img1(:, 2, :);
img1(:, end, :) = img1(:, end-1, :);

img1_ext = reshape(img1(:)*[9/16 3/16 1/16], Nx+2, Ny+2, Nimg, 3);

img2 = zeros(Nx*2, Ny*2, Nimg);
img2(1:2:end, 1:2:end, :) = img1_ext(1:end-2, 1:end-2, :, 3) + img1_ext(2:end-1, 1:end-2, :, 2) + ...
                            img1_ext(1:end-2, 2:end-1, :, 2) + img1_ext(2:end-1, 2:end-1, :, 1);
img2(2:2:end, 1:2:end, :) = img1_ext(2:end-1, 1:end-2, :, 2) + img1_ext(3:end, 1:end-2, :, 3) + ...
                            img1_ext(2:end-1, 2:end-1, :, 1) + img1_ext(3:end, 2:end-1, :, 2);
img2(1:2:end, 2:2:end, :) = img1_ext(1:end-2, 2:end-1, :, 2) + img1_ext(2:end-1, 2:end-1, :, 1) + ...
                            img1_ext(1:end-2, 3:end, :, 3) + img1_ext(2:end-1, 3:end, :, 2);
img2(2:2:end, 2:2:end, :) = img1_ext(2:end-1, 2:end-1, :, 1) + img1_ext(3:end, 2:end-1, :, 2) + ...
                            img1_ext(2:end-1, 3:end, :, 2) + img1_ext(3:end, 3:end, :, 3);
                        
end
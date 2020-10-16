% load img.mat
cut1 = 1150;

img1 = GetBoneEdge_test1(img, cut1);

FOV = 225;
pixelsize = FOV/512;
edgekernel = 5;
edgescale = 2;
img2 = GetBoneEdge_test2(img, cut1, pixelsize, edgekernel, edgescale);
 
% h = 1.0;
% [Nx, Ny] = size(img);
% 
% D1 = (img>cut1).*2.0-1.0;
% D1f = fft2(D1);
% 
% [X, Y] = ndgrid([0:Nx/2 -Nx/2+1:-1]./Nx, [0:Ny/2 -Ny/2+1:-1]./Ny);
% sigma = h/10;
% f1 = exp(-(X.^2+Y.^2)./sigma.^2);
% 
% D2 = ifft2(D1f.*f1);
% D3 = abs(D2(2:end-1,3:end) - D2(2:end-1, 1:end-2)) + abs(D2(3:end, 2:end-1) - D2(1:end-2, 2:end-1))./h;
% 
% D4 = zeros(Nx, Ny);
% D4(2:end-1,2:end-1) = D2(2:Nx-1, 2:Ny-1).*D3;
% D4(D4<0) = 0;
% 
% scale2 = 2.0;
% img2inh = img-cut1;
% img2inh(img2inh<0) = 0;
% img2 = img + img2inh.*D4.*scale2;
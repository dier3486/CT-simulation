function Imageout = GetBoneEdge_test2(ImgIn, minValue, pixelsize, edgekernel, edgescale)
% suggest pixelsize=FOV/512, edgekernel=5-10, edgescale=2

[Nx, Ny] = size(ImgIn);

D1 = (ImgIn>minValue).*2.0-1.0;
D1 = fft2(D1);

[X, Y] = ndgrid([0:Nx/2 -Nx/2+1:-1]./Nx, [0:Ny/2 -Ny/2+1:-1]./Ny);
sigma = pixelsize/edgekernel;
f1 = exp(-(X.^2+Y.^2)./sigma.^2);

D1 = ifft2(D1.*f1);
D2 = abs(D1(2:end-1,3:end) - D1(2:end-1, 1:end-2)) + abs(D1(3:end, 2:end-1) - D1(1:end-2, 2:end-1))./pixelsize;

D3 = zeros(Nx, Ny);
D3(2:end-1,2:end-1) = D1(2:Nx-1, 2:Ny-1).*D2;
D3(D3<0) = 0;

img2inh = ImgIn - minValue;
img2inh(img2inh<0) = 0;
Imageout = ImgIn + img2inh.*D3.*edgescale;
         
end

function ImgOut = bonecorrEdgeEnhance3D(ImgIn, minValue, pixelsize, edgekernel, edgescale)
% enhance the bone

[Ny, Nx, Nz] = size(ImgIn);
if size(pixelsize(:),1) == 1
    pixelsize = repmat(pixelsize, 3, 1);
end

[X, Y, Z] = meshgrid([0:Nx/2 -Nx/2+1:-1]./Nx, [0:Ny/2 -Ny/2+1:-1]./Ny, [0:Nz/2 -Nz/2+1:-1]./Nz);
sigma = pixelsize./edgekernel;
f1 = exp(-((X./sigma(1)).^2+(Y./sigma(2)).^2+(Z./sigma(3)).^2));

D1 = ifftn(fftn((ImgIn>minValue).*1.0).*f1, 'symmetric');
D2 =  abs(D1(3:end, 2:end-1, 2:end-1) - D1(1:end-2, 2:end-1, 2:end-1))./pixelsize(1) + ...
      abs(D1(2:end-1, 3:end, 2:end-1) - D1(2:end-1, 1:end-2, 2:end-1))./pixelsize(2) + ...
      abs(D1(2:end-1, 2:end-1, 3:end) - D1(2:end-1, 2:end-1, 1:end-2))./pixelsize(3);


% D2 = (abs(D1(2:end-1,3:end, :) - D1(2:end-1, 1:end-2, :)) + abs(D1(3:end, 2:end-1, :) - D1(1:end-2, 2:end-1, :)))./pixelsize;

D3 = ImgIn.*0;
D3(2:end-1,2:end-1, 2:end-1) = D1(2:Ny-1, 2:Nx-1, 2:Nz-1).*D2;
D3(D3<0) = 0;

img2inh = ImgIn - minValue;
img2inh(img2inh<0) = 0;
ImgOut = ImgIn + img2inh.*D3.*edgescale;
end
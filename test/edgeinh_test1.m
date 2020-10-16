figh = figure;
D = image;
D = D.CData;
close(figh);

Nx = 64; Ny = 64;
Nx0 = 128;  Ny0 = 128;
D0 = zeros(Nx0, Ny0);
D0(1:64,1:64) = D;

bcut = 25;
D1 = (D0>bcut).*2.0-1.0;
D1f = fft2(D1);


[X, Y] = ndgrid([0:Nx0/2 -Nx0/2+1:-1], [0:Ny0/2 -Ny0/2+1:-1]);
sigma = 1000.0;
f1 = exp(-(X.^2+Y.^2)./sigma);

D2 = ifft2(D1f.*f1);
D3 = abs(D2(2:end-1,3:end) - D2(2:end-1, 1:end-2)) + abs(D2(3:end, 2:end-1) - D2(1:end-2, 2:end-1));

D4 = zeros(Nx, Ny);
D4(2:end-1,2:end-1) = D2(2:Nx-1, 2:Ny-1).*D3(1:Nx-2, 1:Ny-2);
D4(D4<0) = 0;

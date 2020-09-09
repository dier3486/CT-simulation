% input image
image0 = img1(:,:,8);

wincenter = 1055;
winwidth = 100;

Lb = 1055-50;
Ub = 1055+50;

Niter = 100;
imagesize = size(image0, 1);

u = image0;
u0 = u;

for ii = 1:Niter
    ux = [diff(u, 1, 1); zeros(1, imagesize)];
    uy = [diff(u, 1, 2)  zeros(imagesize, 1)];
    
end
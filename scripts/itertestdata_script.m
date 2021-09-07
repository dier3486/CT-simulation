% paramters
pbm.delta_d = 0.5;
pbm.h = 0.6;
pbm.Npixel = 700;
pbm.midchannel = (pbm.Npixel+1)/2;
pbm.Nview = 500;
pbm.imagesize = 512;
pbm.fillmiss = 0;
pbm.maxR = (pbm.Npixel+1)/2*pbm.delta_d/pbm.h;

viewangle = linspace(0, pi, pbm.Nview+1);
pbm.viewangle = viewangle(1:end-1);
pbm = everything2single(pbm);

% phantom
img0 = single(phantom('Modified Shepp-Logan',pbm.imagesize));
img0 = img0.*5.0e3;
% week the shell
img0(img0>2500) = 2500;
% set the 3 dots to lower density
tmp = img0(390:430, 210:290);
tmp(tmp>1000.01) = 1000+(tmp(tmp>1000.01)-1000).*2e-3.*20.0;
img0(390:430, 210:290) = tmp;
% blur a little
Cblur = 1.5;
img0 = gaussblur(img0, Cblur);

% a dot
% img0(382:385,159:162) = 975;

% [Y, X] = meshgrid(-(pbm.imagesize-1)/2:(pbm.imagesize-1)/2);
% img0 = zeros(pbm.imagesize, pbm.imagesize, 'single');
% img0( sqrt(X.^2+Y.^2).*pbm.h < 100) = 1.0;

gpuDevice;

% FP
P0 = parallelprojinimageGPU(pbm, img0);
% or,  P0 = parallelprojinimage(pbm, img0, '2D linearinterp');

% noise
mu_e = 0.001;
Ndose = 4.0e6;
% mu_e = 0.016;
% Ndose = 0.5e6;
I0 = exp(-gather(P0).*1e-3.*mu_e);
I1 = poissrnd(I0.*Ndose)./Ndose;
P1 = -log(I1)./mu_e.*1e3;
P1(P1>1e6) = 1e6;
P1 = gpuArray(P1);

% BP
img1 = filterbackproj2D(P1, pbm, 'hann');

% % noise intervene
% n = 1024;
% c_intv = 0.9;
% noise_intv = linspace(1, 0, n/2+1);
% noise_intv = [noise_intv fliplr(noise_intv(2:end-1))];
% DI = ifft(fft(I1-I0, n, 2).*(1-c_intv*0.8+noise_intv.*c_intv), n, 2);
% I2 = I0+DI(:, 1:pbm.Nview);
% P2 = -log(I2)./mu_e.*1e3;
% P2(P2>1e6) = 1e6;
% P2 = gpuArray(P2);
% img2 = filterbackproj2D(P2, pbm, 'hann');






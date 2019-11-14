% test script for ideal detector position
% on grid image
% 2D

addpath(genpath('../'));

Nview = 1440;
% Nview = 3;
viewangle = linspace(0, pi*2, Nview+1);
viewangle = viewangle(1:end-1);


parallelbeam.Np = 950;
parallelbeam.midchannel = 476.25;
% parallelbeam.midchannel = 474.75;
parallelbeam.delta_d = 0.55;
parallelbeam.h = 500/512;
parallelbeam.viewangle = viewangle;
% parallelbeam.viewangle = mod(viewangle - pi/2, pi*2);
parallelbeam.N = 512;

d = parallelbeam.h/parallelbeam.delta_d;
% M = backproj2D_spM1(parallelbeam.Np, parallelbeam.viewangle, parallelbeam.midchannel, d, parallelbeam.N);

% read Cimage
fid = fopen('D:\matlab\data\simulation\20180828\td\obj.bin');
Cimage = fread(fid,inf,'single=>single');
fclose(fid);
Cimage=reshape(Cimage, 512, 512);

% read filter
fid = fopen('D:\matlab\ct\kernel\BodySoft.bin.res');
myfilter = fread(fid, inf, 'single=>single');
fclose(fid);

% read projection
fid = fopen('D:\matlab\data\simulation\20180828\td\pro.bin');
pro1 = fread(fid,inf,'single=>single');
fclose(fid);
pro1=reshape(pro1, 950, []);

% FP1
tic;
D1 = parallelprojinimage(parallelbeam, Cimage, '2D linearinterp');
toc

% BP1
tic;
[B1,H] = filterbackproj2D(D1, parallelbeam, myfilter);
toc;
% B1 = B1.*2;

% Niter = 10;
% B0 = Cimage;
% alpha = 1.0;
% for ii = 1:Niter
%     tic
%     Dii = parallelprojinimage(parallelbeam, B0, '2D linearinterp');
%     [Bii, H] = filterbackproj2D(Dii, parallelbeam, myfilter);
%     Bii(isnan(Bii)) = 0;
% %     Bii = Bii.*2;
%     Rii = Cimage - Bii;
%     B0 = B0 + Rii.*alpha;
%     toc
% end

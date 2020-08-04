% images
img0 = loaddcmimage('F:\data-Dier.Z\imgdata\701_Original');
% input paramters
imgcenter = [0, 2.1091623];
FOV = 250;
hx_ISO = 0.6110;
imagesize = 512;
Nslice = 32;
pitch = 1/2;

% T
img0T = permute(img0, [2 1 3]);

% anti-ring paramters
ringfilter1 = [-1/2 1 -1/2];
% ringfilter1 = [-1/4 -1/4 1 -1/4 -1/4];
% ringfilter1 = [-1/5 -1/5 4/5 -1/5 -1/5];
Nrow = 16;
% I know Nrow = Nslice*pitch;
rotdirect = 1;
Lb = 1050-150;
Ub = 1050+150;
Ntheta = 192;
d_radius = hx_ISO*imagesize/FOV;
Zsample = 2;
% I know Zsample = 1/pitch;
restcut1 = 0.3;

imgfix1 = antiringonhelical(img0T, imgcenter, Nrow, rotdirect, Lb, Ub, Ntheta, d_radius, Zsample, true, restcut1, ringfilter1);

% fixed image1
img1T = img0T - imgfix1;
img1 = permute(img1T, [2 1 3]);
% try to compare img1 with img0 ...

% merge to thick image
Nmerge = 8;
% anti-ring again on thick image
img1Tm = squeeze(mean(reshape(img1T, imagesize, imagesize, Nmerge, []), 3));
ringfilter2 = [-1/3 2/3 -1/3];
% ringfilter2 = [-1/5 -1/5 4/5 -1/5 -1/5];
% ringfilter2 = [-1 -1 -1 -1 8 -1 -1 -1 -1]./9;
restcut2 = 0.1;
imgfix2 = antiringonimage(img1Tm, imgcenter, Lb, Ub, Ntheta, d_radius, true, restcut2, ringfilter2);
% fixed image2
img2Tm = img1Tm - imgfix2;
img2m = permute(img2Tm, [2 1 3]);

% try to compare img2m with img0m and img1m ...
img0m = squeeze(mean(reshape(img0, imagesize, imagesize, Nmerge, []), 3));
img1m = squeeze(mean(reshape(img1, imagesize, imagesize, Nmerge, []), 3));
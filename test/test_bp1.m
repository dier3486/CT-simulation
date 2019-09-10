% test script for ideal detector position
% on grid image
% 2D

addpath(genpath('../'));

img = zeros(100,100);
img(25:75, 25:75) = 1;
Nview = 200;
theta = linspace(0, pi*2, Nview+1);
theta = theta(1:end-1);
theta_180 = theta.*(180/pi);
[R0,xp0] = radon(img, theta_180);

prm.Np = 201;
prm.midchannel = 101;
prm.delta_d = 0.7;
prm.h = 1;
prm.viewangle = theta;
prm.imagesize = 100;

P1 = parallelprojinimage(prm, img, '2D linearinterp');

d = prm.delta_d/prm.h;
B1 = iradon(P1, theta_180, d, 100);
[B2,H1] = filterbackproj2D(P1, prm);
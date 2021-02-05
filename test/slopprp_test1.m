
Nslice = 128;
hz = 0.6222;

Nviewprot = 1152;
SID = 570;
theta_tilt = 0.2;
imagesize = 512;
delta_xy = 0.55;
FOV = 300;

thetaview = (0:Nviewprot-1).*(pi*2/Nviewprot);
[X, Y] = ndgrid((-imagesize+1)/2 : (imagesize-1)/2);

ZgridS = -Nslice/2+1:Nslice/2;
ZgridA = reshape(-Nslice+1:Nslice, 1, 1, []);

r = sqrt(X.^2+Y.^2);
theta_xy = atan2(Y, X);

n1 = 17;

% Nv = 300;
Nv = 60;
Vs = linspace(0, pi*2, Nv+1);
% Nr = 256;
Nr = 16;
Rs = linspace(0, FOV/2, Nr)';

eta1 = Rs*sin(Vs)./SID;
zeta1 = Rs*cos(Vs)./SID;
t1 = axialconehomeomorph(eta1, zeta1, Nslice, theta_tilt);

Nez = 33;
ezgrid = linspace(-FOV/2/SID, FOV/2/SID, Nez);
[eta2, zeta2] = ndgrid(ezgrid);
t2 = axialconehomeomorph(eta2, zeta2, Nslice, theta_tilt);



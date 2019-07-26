% test script for ideal detector position

addpath(genpath('../'));

% system parameters (sample)
% focalpos_ideal = [0, -550, 0];
SSD = 550;
SDD = 1000;
hx_ISO = 0.55;
hz_ISO = 0.55;
Npixel = 950;
mid_U = 475.25;     % satrt from 1
Nslice = 16;

% detectors position
alpha_1 = asin(hx_ISO/SSD);
alpha_pixel = ((1:Npixel)-mid_U).*alpha_1 + pi/2;
x_pos = cos(alpha_pixel).*SDD;
y_pos = sin(alpha_pixel).*SDD - SSD;
z_pos = ((1:Nslice)-Nslice/2-1/2).*(hz_ISO*SDD/SSD);
XYZ = [repmat(x_pos(:), Nslice, 1), repmat(y_pos(:), Nslice, 1), repelem(z_pos(:), Npixel, 1)];



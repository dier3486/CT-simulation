% test script for ideal detector position
% on grid image

addpath(genpath('../'));

% image
Cimage0 = reshape(1:256, 16, 16);
% Cimage = ones(16, 16);
h = 32.0;

% parallel beams
Np = 20;
delta_d = 30;
mid_chn = 10.25;
% view angles
viewangle = 1.1+pi/2;

% size
[Nx, Ny] = size(Cimage0);
Xgrid = (-Nx/2:Nx/2).*h;
Ygrid = (-Ny/2:Ny/2).*h; 

Nc = Nx*Ny;
Cimage = [Cimage0(:); 0];
d = ((1:Np)-mid_chn).*delta_d;
d = d(:);
d_h = d./h;

% suppose pi/4 < viewangle < pi*3/4



dt_x = tan(pi/2-viewangle);

t_x = repmat((-(Ny-1)/2:(Ny-1)/2).*dt_x, Np, 1) - repmat(d_h.*csc(viewangle), 1, Ny) + 1/2;

index_x1 = floor(t_x);
inter_alpha = t_x - index_x1;

index_x1 = index_x1 + Nx/2;
index_x2 = index_x1 + 1;
index_x1(index_x1<=0 | index_x1>Nx) = nan;
index_x2(index_x2<=0 | index_x2>Nx) = nan;
index_x1 = index_x1 + repmat((0:Ny-1).*Nx, Np, 1);
index_x2 = index_x2 + repmat((0:Ny-1).*Nx, Np, 1);
index_x1(isnan(index_x1)) = Nc+1;
index_x2(isnan(index_x2)) = Nc+1;

dt = abs(h./sin(viewangle));
P1 = sum(Cimage(index_x1).*(1-inter_alpha) + Cimage(index_x2).*inter_alpha, 2).*dt;

[inter_a, index_1, index_2, cs_view] = parallellinearinterp2D(Nx, Ny, d_h, viewangle, Nc);
P1b = sum(Cimage(index_1).*(1-inter_a) + Cimage(index_2).*inter_a, 2).*(abs(cs_view)*h);

theta = repmat(viewangle, Np, 1);
[Vdt, Vindex] = linesinimage2D(theta, d, [], 0, Xgrid, Ygrid);
P2 = sum(Vdt.*Cimage(Vindex), 2);

parallelbeam.Np = Np;
parallelbeam.delta_d = delta_d;
parallelbeam.midchannel = mid_chn;
parallelbeam.h = h;

P3 = parallelprojinimage(parallelbeam, Cimage0, viewangle);


% bone projection test1


gpuDevice;

% FP
Nview = 400;
deltaview = pi/Nview; 
viewangle = linspace(0, pi, Nview+1);
viewangle = single(viewangle(1:Nview));

maxFOV = 500;
h_img = prmflow.recon.FOV/prmflow.recon.imagesize;
dp = h_img*1.8;

Np = floor(min(prmflow.recon.FOV*sqrt(2), maxFOV)/dp/2)*2+1;
d_h = gpuArray(single((-(Np-1)/2:(Np-1)/2)'.*(dp/h_img)));
Cimage = gpuArray(dataflow.image);
zetaangle = gpuArray(viewangle(1:Nview/2) + pi/2);
imgcenter = [0 0];

D = parallelcrossprojectimg2D(Cimage, d_h, zetaangle, imgcenter);




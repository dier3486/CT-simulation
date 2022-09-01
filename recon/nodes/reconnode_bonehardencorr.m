function [dataflow, prmflow, status] = reconnode_bonehardencorr(dataflow, prmflow, status)
% recon node, bone-beamharden correction
% [dataflow, prmflow, status] = reconnode_bonehardencorr(dataflow, prmflow, status); 

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% parameters to use in prmflow
recon = prmflow.recon;
imagesize = recon.imagesize;
voxelsize = recon.voxelsize;
delta_d = recon.delta_d;
hond = voxelsize/delta_d;
Npixel = recon.Npixel;
midchannel = recon.midchannel;
Nviewprot = recon.Nviewprot;
reconcenter = recon.center;
delta_view = recon.delta_view;
startviewangle = recon.startviewangle;
Nimage = recon.Nimage;
effFOV = recon.effFOV;
effNp = floor(effFOV/delta_d) + 1;

% calibration table
bonecorr = prmflow.corrtable.(status.nodename);
bonecurve = reshape(bonecorr.bonecurve, bonecorr.bonecurvelength, []);

nodeprm = prmflow.pipe.(status.nodename);
% enhance bone edge
if isfield(nodeprm, 'edgeenhance')
    edgeenhance = nodeprm.edgeenhance;
else
    edgeenhance = true;
end
if isfield(nodeprm, 'edgekernel')
    edgekernel = nodeprm.edgekernel;
else
    edgekernel = 2.0;
end
if isfield(nodeprm, 'edgescale')
    edgescale = nodeprm.edgescale;
else
    edgescale = 1.0;
end
if edgeenhance
    minvalue = min(bonecurve(:,1)) - 20;
    pixelsize = prmflow.recon.FOV/prmflow.recon.imagesize;
    ImageBE = bonecorrEdgeEnhance(dataflow.image, minvalue, pixelsize, edgekernel, edgescale);
else
    ImageBE = dataflow.image;
end

% use bone curve to adjust image
BoneImage = GetBoneImg(ImageBE, bonecurve);

% general perpare

% sub view
if isfield(nodeprm, 'subview')
    subview = nodeprm.subview;
else
    subview = 1;
end
viewangle = mod((0:Nviewprot/2-1).*delta_view + startviewangle(1) + pi/2, pi*2);
viewangle = viewangle(1:subview:end);
Nview = length(viewangle);
eta_C = reconcenter(1).*sin(viewangle) - reconcenter(2).*cos(viewangle);
indexstart = floor(midchannel + (-effFOV/2 + eta_C)./delta_d);
indexstart(indexstart<1) = 1;
ctrIdx = midchannel+eta_C./delta_d+1-indexstart;
channelpos = ((1:Npixel)'-midchannel).*delta_d;
maxR = effFOV/2/delta_d;

% BBH table
HCscale = 1000;
Dscale = gpuArray(1/HCscale/bonecorr.curvescale(1));
bonescale = gpuArray(1/HCscale/bonecorr.curvescale(2));
bonecorr.order = reshape(bonecorr.order, 1, []);
curvematrix = gpuArray(reshape(bonecorr.curvematrix, bonecorr.order));
efffilter = interp1( bonecorr.beamposition, bonecorr.effbeamfilter, channelpos, 'linear', 'extrap');
efffilter = gpuArray(efffilter./bonecorr.curvescale(3));
mubonmuw = gpuArray(single(bonecorr.refrencebonemu/bonecorr.refrencemu));
Nw = 400; Nb = 200; Nf = 20;
gridW = gpuArray(single(linspace(0, 0.5, Nw)));
gridB = gpuArray(single(linspace(0, 2.0, Nb)));
gridF = gpuArray(single(linspace(0.1, 1.0, Nf)));
[gBB, gWW, gFF] = meshgrid(gridB, gridW, gridF);
BBHmatrix = polyval3dm(curvematrix, gWW, gBB, gFF).*mubonmuw-1;

% gridW = gpuArray(gridW);
% gridB = gpuArray(gridB);
% gridF = gpuArray(gridF);
% BBHmatrix = gpuArray(BBHmatrix);

% Filter
if isfield(nodeprm, 'Filter')
    filter = gpuArray(loadfilter(nodeprm.Filter, prmflow.recon.Npixel, prmflow.recon.delta_d));
else
    filter = gpuArray(prmflow.recon.filter);
    if isfield(recon, 'upsampling') && recon.upsampling
        warning('The filter in boneharden correction coulde be wrong!');
    end
end
Hlen = length(filter);

% BP prepare
[x,y,z,Sxy,costheta,sintheta] = backproj2Dprepare(imagesize, Nimage, hond, [0 0], maxR, viewangle, 'single');
Nxy = size(x, 1);
ctrIdx = gpuArray(cast(ctrIdx, 'single'));
image_fix = zeros(Nxy, Nimage, 'single', 'gpuArray');

% FP prepare
Nx = gpuArray(single(imagesize));
Ny = Nx;
d_h = gpuArray(single(channelpos./voxelsize));
image0 = gpuArray(dataflow.image);
BoneImage = gpuArray(BoneImage);
viewangleGPU = gpuArray(viewangle);
imgindex = gpuArray(repmat(reshape(single(1:Nimage), 1, 1, []), effNp, Nx));
reconcenter_h = gpuArray(single(reconcenter./voxelsize));

% P0 = zeros(effNp, Nimage, Nview, 'single');

% tic; PImageBE = FP(image0, Nx, Ny, Np, Nviewprot, h_img, d_h); toc;

% 1;
% tic;
for iview = 1:Nview
    dh_iview = d_h(indexstart(iview):indexstart(iview)+effNp-1);
    effF_ivew = repmat(efffilter(indexstart(iview):indexstart(iview)+effNp-1), 1, Nimage);
    % warn: the error of the efffilter by the images center movement due to
    % the gantry tilt is ignored.
    [interpX, interpY, cs_view] = parallellinearinterp2D2(Nx, Ny, dh_iview, viewangleGPU(iview), reconcenter_h);
    % I know Nx = Ny, if not these interp3 will catch a bug in size of imgindex.
    interpY_rep = repmat(interpY, 1, 1, Nimage);
    interpX_rep = repmat(interpX, 1, 1, Nimage);
    D0 = squeeze(sum(interp3(image0, interpY_rep, interpX_rep, imgindex, 'linear', 0), 2)).*(abs(cs_view)*voxelsize);
    DB = squeeze(sum(interp3(BoneImage, interpY_rep, interpX_rep, imgindex, 'linear', 0), 2)).*(abs(cs_view)*voxelsize);
    D0 = D0.*(D0>0);
    DB = DB.*(DB>0);
    % BBH
    % I thought it is too slow
%     Dbcurve = polyval3dm(curvematrix, D0.*Dscale, DB.*bonescale, efffilter(indexstart(iview):indexstart(iview)+effNp-1));
%     Dfix0 = (Dbcurve.*mubonmuw-1).*DB;
    Dfix = interp3(gBB, gWW, gFF, BBHmatrix, DB.*bonescale, D0.*Dscale, effF_ivew).*DB;
    
    % debug
    if any(isnan(Dfix(:)))
        1;
    end
    % filter
    Df = ifft(fft(Dfix, Hlen).*filter, 'symmetric');
    
%     % debug
%     Df = ifft(fft(D0, Hlen).*filter, 'symmetric');

    % BP
    Eta = (-x.*sintheta(iview) + y.*costheta(iview)) + ctrIdx(iview);
    image_fix = image_fix + interp2(Df, z, Eta, 'linear', 0);
    
%     P0(:,:,iview) = gather(D0);
end
% toc;

% image diff
image_diff = zeros(imagesize^2, Nimage, 'single');
image_diff(Sxy, :) = gather(image_fix).*(pi/Nview/2);
image_diff = reshape(image_diff, imagesize, imagesize, Nimage);

% add diff to original
dataflow.image = dataflow.image + image_diff;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function ImgOut = GetBoneImg(ImgIn, BoneCurve)
minValue = min(BoneCurve(:,1));
maxValue = max(BoneCurve(:,1));
ImgIn(ImgIn < minValue) = 0;
ImgIn(ImgIn > maxValue) = maxValue;
idx = find(ImgIn > 0);
ImgIn(idx)=interp1(BoneCurve(:,1), BoneCurve(:,2), ImgIn(idx));
ImgOut = ImgIn;
end


function [x,y,z,Sxy,costheta,sintheta] = backproj2Dprepare(N, Nslice, hond, centerond, maxR, theta, pclass)

% Define the x & y axes for the reconstructed image
[x, y] = ndgrid(-(N-1)/2 : (N-1)/2);
x = x(:).*hond - centerond(:, 1)';
y = y(:).*hond - centerond(:, 2)';
if size(x, 2)==1 && Nslice>1
    x = repmat(x, 1, Nslice);
    y = repmat(y, 1, Nslice);
end
Sxy = any(x.^2 + y.^2 <= maxR.^2, 2);

x = x(Sxy, :);
y = y(Sxy, :);
Nxy = sum(Sxy);
% z (slice)
z = repmat(1:Nslice, Nxy, 1);

% Generate trignometric tables
costheta = cos(theta(:)');
sintheta = sin(theta(:)');

% to gpu
x = gpuArray(cast(x, pclass));
y = gpuArray(cast(y, pclass));
z = gpuArray(cast(z, pclass));
costheta = gpuArray(cast(costheta, pclass));
sintheta = gpuArray(cast(sintheta, pclass));
% ctrIdx = gpuArray(cast(ctrIdx, pclass));
end
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

% GPU?
if isfield(status, 'GPUinfo') && ~isempty(status.GPUinfo)
    GPUonoff = true;
else
    GPUonoff = false;
end
% parameters to use in prmflow
recon = prmflow.recon;
imagesize = recon.imagesize;
voxelsize = recon.voxelsize;
delta_d = recon.delta_d;
Npixel = recon.Npixel;
midchannel = recon.midchannel;
Nviewprot = recon.Nviewprot;
Nimage = recon.Nimage;
reconcenter_h = recon.center./recon.voxelsize;
FPchannelpos = recon.FPchannelpos;
FPchannelpos_h = FPchannelpos./recon.voxelsize;
effNp = min(recon.effNp, Npixel);

% parameters in pipe
nodeprm = prmflow.pipe.(status.nodename);

% viewangle
% sub view
if isfield(nodeprm, 'subview')
    subview = nodeprm.subview;
else
    subview = 1;
end
viewangle = mod(recon.viewangle(1:subview:Nviewprot/2), pi*2);
indexstart = recon.indexstart(1:subview:Nviewprot/2, 1);
Nview = length(viewangle);

% calibration table
bonecorr = prmflow.corrtable.(status.nodename);
bonecurve = reshape(bonecorr.bonecurve, bonecorr.bonecurvelength, []);
% BBH table
bonecorr = loadboneharden(bonecorr, FPchannelpos);
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
    ImageBE = bonecorrEdgeEnhance(dataflow.image, minvalue, voxelsize, edgekernel, edgescale);
else
    ImageBE = dataflow.image;
end
% use bone curve to adjust image
BoneImage = getBoneImg(ImageBE, bonecurve);

% Filter
if isfield(nodeprm, 'Filter')
    filter = loadfilter(nodeprm.Filter, recon.Npixel, recon.delta_d);
else
    filter = prmflow.recon.filter;
    if isfield(recon, 'upsampling') && recon.upsampling
        warning('The filter in boneharden correction coulde be wrong!');
    end
end

% BP prepare
XY_d = recon.XY.*(recon.SID/delta_d);
Sxy = recon.activeXY;
Nxy = recon.NactiveXY;
Zindex = repmat(single(1:Nimage), Nxy, 1);
image_fix = zeros(Nxy, Nimage, 'single');

% FP prepare
image0 = dataflow.image + BoneImage.*1i;
imgindex = repmat(reshape(single(1:Nimage), 1, 1, []), effNp, max(imagesize));

% to GPU
if GPUonoff
    [XY_d, FPchannelpos_h, indexstart, effNp, midchannel, imgindex, voxelsize, image0, image_fix, ...
     Zindex, bonecorr, filter, viewangle, Nimage, imagesize] = putinGPU...
    (XY_d, FPchannelpos_h, indexstart, effNp, midchannel, imgindex, voxelsize, image0, image_fix, ...
     Zindex, bonecorr, filter, viewangle, Nimage, imagesize);
end

Hlen = length(filter);
sintheta = sin(viewangle);
costheta = cos(viewangle);
D0max = max(bonecorr.Wmesh(:))/bonecorr.Dscale;
DBmax = max(bonecorr.Bmesh(:))/bonecorr.bonescale;
mubonmuw = bonecorr.refrencebonemu/bonecorr.refrencemu;

% tic;
for iview = 1:Nview
    dh_iview = FPchannelpos_h(indexstart(iview):indexstart(iview)+effNp-1);
    effF_ivew = repmat(bonecorr.efffilter(indexstart(iview):indexstart(iview)+effNp-1), 1, Nimage);
    % warn: the error of the efffilter by the images center movement due to
    % the gantry tilt is ignored.
    [interpX, interpY, cs_view] = parallellinearinterp2D2(imagesize(1), imagesize(2), dh_iview, ...
        viewangle(iview), reconcenter_h);
    interpX_rep = repmat(interpX, 1, 1, Nimage);
    interpY_rep = repmat(interpY, 1, 1, Nimage);
    Nsmp = size(interpX, 2);
    D = squeeze(sum(interp3(image0, interpX_rep, interpY_rep, imgindex(:, 1:Nsmp, :), 'linear', 0), 2)).* ...
         (abs(cs_view)*voxelsize);
    D0 = real(D);
    DB = imag(D);
    D0(D0>D0max) = D0max;
    DB(DB>DBmax) = DBmax;
    % BBH
%     Dfix = (polyval3dm(bonecorr.curvematrix, D0.*bonecorr.Dscale, DB.*bonecorr.bonescale, effF_ivew).*mubonmuw-1).*DB;
    % I thought the interp3 could be faster than polyval.
    Dfix = interp3(bonecorr.Bmesh, bonecorr.Wmesh, bonecorr.Fmesh, bonecorr.BBHmatrix, ...
                   DB.*bonecorr.bonescale, D0.*bonecorr.Dscale, effF_ivew, 'linear', 0).*DB;
    
    % debug
    if any(isnan(Dfix(:)))
        1;
    end
    % filter
    Df = ifft(fft(Dfix, Hlen).*filter, 'symmetric');
    
    % BP
    Eta = -XY_d(:,1).*sintheta(iview) + XY_d(:,2).*costheta(iview);
    t_chn = repmat(Eta + midchannel + 1 - indexstart(iview), 1, Nimage);
    image_fix = image_fix + interp2(Df, Zindex, t_chn, 'linear', 0);
end
% toc;

% image diff
image_diff = zeros(imagesize(1)*imagesize(2), Nimage, 'single');
image_diff(Sxy, :) = gather(image_fix).*(pi/Nview/2);
image_diff = reshape(image_diff, imagesize(2), imagesize(1), Nimage);

% add diff to original
dataflow.image = dataflow.image + image_diff;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end

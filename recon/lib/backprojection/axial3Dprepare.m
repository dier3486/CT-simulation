function recon = axial3Dprepare(recon, BPprm)
% more prepare works for 3D Axial

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

% Nimage and images center
recon.Nimage = recon.Nslice * recon.Nshot;
[recon.imagecenter, recon.reconcenter_2DBP] = imagescenterintilt(recon.center, recon);

% recon range
reconD = sqrt(sum((recon.FOV/2+abs(recon.center)).^2))*2;
% Neighb and Nextslice
if isfield(BPprm, 'Neighb')
    Neighb = BPprm.Neighb;
else
    Rfov = min(sqrt(sum(recon.center.^2)) + recon.effFOV/2, recon.maxFOV/2);
    Rfov = min(Rfov, recon.maxFOV/2);
    Neighb = floor((recon.Nslice*recon.delta_z/2 - (sqrt(recon.SID^2-Rfov^2) - Rfov)/recon.SID*(recon.Nslice-1) ...
             /2*recon.delta_z)/recon.imageincrement) + 2;
end
Neighb = min(Neighb, recon.Nslice/2);
recon.Neighb = Neighb;
recon.Nextslice = recon.Nslice + Neighb*2;

Npixel = recon.Npixel;
% is X-upsampling?
if recon.upsampling
    midchannel_up = recon.midchannel_up;
    Npixel_up = Npixel*2;
else
    midchannel_up = recon.midchannel;
    Npixel_up = Npixel;
end
% Z cross sampling
x1 = max(1, -Npixel_up+midchannel_up*2);
x1 = ceil((x1+1)/2)*2-1;
x2 = min(Npixel_up, midchannel_up*2-1);
x2 = floor(x2/2)*2;
recon.crosschannel = [x1 x2];
recon.Ncrosschn = x2-x1+1;
% Nchannel = prmflow.recon.Ncrosschn;
% Note: the Z cross sampling could be usable while the recon.upsampling is false.

% Z cross interp (in Z rebin)
if isfield(BPprm, 'ZcrossGamma') && ~isempty(BPprm.ZcrossGamma)
    ZcrossGamma = BPprm.ZcrossGamma;
else
    ZcrossGamma = [1.0 2.0];
end
recon.Zcrossinterp = Zcrossinterp(recon, ZcrossGamma);

% Zinterp table (in BP)
if isfield(BPprm, 'Zinterptablesize')
    tablesize = BPprm.Zinterptablesize;
else
    tablesize = 512;
end
% coneflag is 2
coneflag = 2;
% set coneflag = 1 to employ old version of cross BP
recon.Zinterp = ZetaEta2TzTable(tablesize, recon.maxFOV, reconD, recon.SID, recon.Nslice, recon.gantrytilt, coneflag);

% Z upsampling
recon.Zupsamp = Zupsamplingprepare(recon, BPprm, coneflag);

% forward projection 
SID_h = recon.SID/recon.voxelsize;
FPchannelpos_h = recon.forward.FPchannelpos./recon.voxelsize;
recon.forward.FPfocal = [[-sqrt(SID_h^2 - FPchannelpos_h.^2); sqrt(SID_h^2 - FPchannelpos_h.^2)] ...
    [FPchannelpos_h; -FPchannelpos_h]];
recon.forward.FPvector = [zeros(recon.Npixel*2, 1)  [FPchannelpos_h; -FPchannelpos_h]] - recon.forward.FPfocal;
% to use the parallellinearinterp3D2.m

end

function [imagecenter, reconcenter_2DBP] = imagescenterintilt(reconcenter, recon)
% Cin is recon center; Cout is the rotation center on images
% 'small-step' style for tilt shots 

imagecenter = repmat(-reconcenter(:)', recon.Nimage, 1);
% Y shfit
% Yshift = -(recon.imageincrement*tan(recon.gantrytilt)).*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
Yshift = -(recon.imageincrement*sin(recon.gantrytilt)).*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
if recon.couchdirection > 0
    Yshift = fliplr(Yshift);
end
Yshift = repmat(Yshift(:), recon.Nshot, 1);
reconcenter_2DBP = imagecenter;
reconcenter_2DBP(:, 2) = reconcenter_2DBP(:, 2) + Yshift;
% Z shift
% Zshift = (recon.imageincrement*sec(recon.gantrytilt)).*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
Zshift = recon.imageincrement.*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
if recon.couchdirection > 0
    Zshift = fliplr(Zshift);
end
Zshift = Zshift(:) - (0:recon.Nshot-1).*(recon.imageincrement.*recon.Nslice).*recon.couchdirection;
Zshift = Zshift(:) - recon.startcouch;
imagecenter = [imagecenter Zshift];
reconcenter_2DBP = [reconcenter_2DBP Zshift];

end

function interpstr = Zcrossinterp(recon, Gamma)

Nslice = recon.Nslice;
% is X-upsampling?
if recon.upsampling
    midchannel_up = recon.midchannel_up;
    delta_d_up = recon.delta_d_up;
else
    midchannel_up = recon.midchannel;
    delta_d_up = recon.delta_d;
end
% Note: the Z cross sampling could be usable while the recon.upsampling is false.

SID = recon.SID;
gantrytilt = recon.gantrytilt;
Ncrosschn = recon.Ncrosschn;
Rslice2img = recon.imageincrement/recon.delta_z;
% Zdistort = 0;   % debug

eta = ((recon.crosschannel(1):recon.crosschannel(2))-midchannel_up).*(delta_d_up/SID);
Z_eta = sqrt(1-eta.^2)./cos(gantrytilt);
Z_eta(abs(eta)<sin(gantrytilt)) = 1;
gaponZ = Nslice./Z_eta.*Rslice2img - Nslice + 1;

interpstr.t = ((1:Nslice/2)-1/2)./Z_eta(:).*Rslice2img + 3/2;
interpstr.gap = (interpstr.t-Nslice/2-1)./gaponZ(:);
s = interpstr.t>Nslice/2+1;
interpstr.t(s) = interpstr.gap(s) + Nslice/2+1;
interpstr.t = [ones(Ncrosschn, 1).*2  interpstr.t+1 Nslice+4-fliplr(interpstr.t)  ones(Ncrosschn, 1).*(Nslice+3)];
% interpstr.t = [ones(Ncrosschn, 1).*2  interpstr.t+1+Zdistort Nslice+4-fliplr(interpstr.t)+Zdistort  ones(Ncrosschn, 1).*(Nslice+3)];

[interpstr.t_odd, interpstr.t_even, interpstr.gamma] = omiga4interp(interpstr.t, Gamma);

end
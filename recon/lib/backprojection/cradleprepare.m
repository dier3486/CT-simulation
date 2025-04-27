function recon = cradleprepare(recon, nodeprm)
% more prepare works for Cradle

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

% Z upsampling
recon.Zupsamp = Zupsamplingprepare(recon, nodeprm);

% Cone weight
if isfield(nodeprm, 'ConeWeightScale')
    recon.ConeWeightScale = nodeprm.ConeWeightScale;
else
    recon.ConeWeightScale = 1.0;
end

% imagesnumber per pitch
recon.imagesperpitch = recon.pitchlength/recon.imageincrement;

% the governing of the images related with the views
% effRadius is the reffective BP radius (to ISO) on SID
recon.effRadius = min(sqrt(sum(recon.center.^2)) + recon.effFOV/2, recon.maxFOV/2)/recon.SID;
if isfield(nodeprm, 'imagenumber') && ~isnan(nodeprm.imagenumber)
    recon.imagenumber = nodeprm.imagenumber;
else
    recon.imagenumber = [];
end
if isfield(nodeprm, 'startview2image')
    recon.startview2image = nodeprm.startview2image / recon.imageincrement;
else
    recon.startview2image = [];
end
% startview2image is the Z distance from the focal position of first view to the first image center along the images increment
% direction, whose unit is image number (/imageincrement).
if isfield(nodeprm, 'startview2bypass')
    recon.startview2bypass = nodeprm.startview2bypass / recon.imageincrement;
else
    recon.startview2bypass = [];
end
% startview2bypass is the Z distance from the focal position of first view to one if an image shall be located.

% viewextra in all view and pi-line recon condition
[recon.viewextra_full, recon.viewextra_pi] = helicalviewexta(recon.Nslice, recon.pitch, recon.effRadius, recon.Nviewprot);

% the governing of the images related with the views
recon.Nviewskip = floor(asin(recon.effRadius)*recon.Nviewprot/pi/2);

% imagerely
recon.imagerely = ceil(recon.viewextra_full / recon.Nviewprot * recon.imagesperpitch);

% update Nview (Nviewskip has been read) 
recon.viewnumber = recon.viewnumber - recon.Nviewskip;
recon.Nview = recon.Nview - recon.Nviewskip;

% forward projection 
SID_h = recon.SID/recon.voxelsize;
FPchannelpos_h = recon.forward.FPchannelpos./recon.voxelsize;
Phi = asin(FPchannelpos_h./SID_h)./(pi*2);
Zgrid_fp = -(recon.Nslice-1)/2 : (recon.Nslice-1)/2;
recon.forward.ZFPa = (sqrt(SID_h^2 - FPchannelpos_h.^2)./SID_h)*Zgrid_fp - Phi.*recon.imagesperpitch;
recon.forward.ZFPb = Zgrid_fp./SID_h;
% to use the parallellinearinterp3D3.m

end
function [dataflow, prmflow, status] = reconnode_Antiring(dataflow, prmflow, status)
% recon node, anti ring on image
% [dataflow, prmflow, status] = reconnode_Antiring(dataflow, prmflow, status);

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

% parameters from pipe
ARprm = prmflow.pipe.(status.nodename);
if isfield(ARprm, 'alpha')
    alpha = ARprm.alpha;
else
    alpha = 1.0;
end
if isfield(ARprm, 'Ntheta')
    Ntheta = ARprm.Ntheta;
else
    Ntheta = 192;
end

% parameters from recon
if isfield(prmflow.recon, 'voxelsize')
    voxelsize = prmflow.recon.voxelsize;
else
    voxelsize = prmflow.recon.FOV/min(prmflow.recon.imagesize);
end
imagecenter = prmflow.recon.imagecenter;
if isfield(prmflow.recon, 'windowcenter') && isfield(prmflow.recon, 'windowwidth')
    Lb =  prmflow.recon.windowcenter - prmflow.recon.windowwidth/2 + 1000;
    Ub =  prmflow.recon.windowcenter + prmflow.recon.windowwidth/2 + 1000;
else
    Lb = 950;
    Ub = 1150;
end
d_radius = prmflow.recon.delta_d/voxelsize;

% cneter on image space
centerfix = imagecenter(:, 1:2)./voxelsize;
% anti ring on image space
imgfix = antiringonimage(dataflow.image, centerfix, Lb, Ub, Ntheta, d_radius);
% add to image
dataflow.image = dataflow.image - imgfix.*alpha;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
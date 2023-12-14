function [prmflow, status] = reconnode_fanangleprepare(prmflow, status)
% prepare node, fananglecorr prepare
% [prmflow, status] = reconnode_fanangleprepare(prmflow, status);

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
Npixel = prmflow.raw.Npixel;
Nslice = prmflow.raw.Nslice;
Nfocal = prmflow.raw.Nfocal;
focalspot = prmflow.raw.focalspot;

% detector
detector = prmflow.system.detector;
if length(detector.mid_U) > 1
    % >= v1.11
    mid_U = detector.mid_U(focalspot);
else
    % v1.0
    mid_U = detector.mid_U;
end
SID = detector.SID;
SDD = detector.SDD;
hx = detector.hx_ISO;
if isfield(detector, 'deltafan') && detector.deltafan ~= 0
    % >= v1.11
    deltafan = detector.deltafan;
else
    % v1.0
    deltafan = atan(mod(mid_U(1), 1)*hx/SID) + atan((1-mod(mid_U(1), 1))*hx/SID);
end

% focalposition
if isfield(detector, 'focalposition')
    focalposition = detector.focalposition(focalspot, :);
else
    focalposition = prmflow.system.focalposition(focalspot, :);
end

[fanangles, focalangle] = detpos2fanangles(detector.position, focalposition);
fanangles = reshape(fanangles, Npixel, Nslice*Nfocal);

% equal angles
if length(detector.mid_U) > 1
    % >= v1.11
    equalfan = ((1:Npixel)' - mid_U(:)').*deltafan + focalangle;
else
    % v1.0
    focal_fix = atan(tan(focalangle-pi/2).*(SID/SDD)) + pi/2;
    equalfan = ((1:Npixel)' - mid_U).*deltafan + focal_fix;
end
equalfan = reshape(repmat(equalfan, Nslice, 1), Npixel, Nslice*Nfocal);

equalfaninterp = zeros(Npixel, Nslice*Nfocal, 'single');
fanindex = single(1:Npixel)';
% interp prepare
for ii = 1:Nslice*Nfocal
    equalfaninterp(:, ii) = interp1(fanangles(:, ii), fanindex, equalfan(:, ii), 'linear', 'extrap');
end
equalfaninterp(equalfaninterp < 1) = 1;
equalfaninterp(equalfaninterp > Npixel) = Npixel;

% to rebin
prmflow.rebin.equalfaninterp = equalfaninterp;
prmflow.rebin.fanangles = equalfan;
prmflow.rebin.focalangle = focalangle;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
function [dataflow, prmflow, status] = reconnode_fananglecorr(dataflow, prmflow, status)
% recon node, interp the rawdata to equal-fanangle
% [dataflow, prmflow, status] = reconnode_fananglecorr(dataflow, prmflow, status);

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
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nfocal = prmflow.recon.Nfocal;
focalspot = prmflow.recon.focalspot;

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

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice*Nfocal, []);

% interp1
for ii = 1:Nslice*Nfocal
    [index_ii, alpha_ii] = interpprepare(fanangles(:, ii), equalfan(:, ii), 'extrap');
    dataflow.rawdata(:, ii, :) = dataflow.rawdata(index_ii(: ,1), ii, :).*alpha_ii(:, 1) + ...
                                 dataflow.rawdata(index_ii(: ,2), ii, :).*alpha_ii(:, 2);
end

% to return
prmflow.recon.fanangles = reshape(equalfan, Npixel*Nslice, Nfocal);
prmflow.recon.focalangle = focalangle;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
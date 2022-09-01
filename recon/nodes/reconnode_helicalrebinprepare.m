function [prmflow, status] = reconnode_helicalrebinprepare(prmflow, status)
% recon node, rebin prepare
% [prmflow, status] = reconnode_helicalrebinprepare(prmflow, status);

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
Nviewprot = prmflow.recon.Nviewprot;
focalspot = prmflow.recon.focalspot;
gantrytilt = prmflow.recon.gantrytilt;
% gantrytilt = prmflow.protocol.gantrytilt*(pi/180);
focalposition = prmflow.system.focalposition(focalspot, :);
% Nfocal = prmflow.recon.Nfocal;
rebinpipe = prmflow.pipe.(status.nodename);

if isfield(rebinpipe, 'viewblock') && ~isempty(rebinpipe.viewblock)
    viewblock = rebinpipe.viewblock;
else
    viewblock = 32;
end

% detector
detector = prmflow.system.detector;

% fan angles & focal angle(s)
if isfield(prmflow.recon, 'fanangles')
    fanangles = prmflow.recon.fanangles;
    focalangle = prmflow.recon.focalangle;
else
    [fanangles, focalangle] = detpos2fanangles(detector.position, focalposition);
%     prmflow.recon.fanangles = reshape(fanangles, prmflow.recon.Npixel, prmflow.recon.Nslice);
    fanangles = reshape(fanangles, prmflow.recon.Npixel, prmflow.recon.Nslice, []);
    prmflow.recon.fanangles = fanangles;
    prmflow.recon.focalangle = focalangle;
end

% rebin prepare
rebin = helicalrebinprepare(detector, fanangles, focalangle, Nviewprot);
rebin.viewblock = viewblock;
% rebin.FOV = prmflow.recon.FOV;

% rebin to prmflow
prmflow.rebin = rebin;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
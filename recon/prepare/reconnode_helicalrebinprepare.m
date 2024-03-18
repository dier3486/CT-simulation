function [dataflow, prmflow, status] = reconnode_helicalrebinprepare(dataflow, prmflow, status)
% recon node, rebin prepare
% [dataflow, prmflow, status] = reconnode_helicalrebinprepare(dataflow, prmflow, status);

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

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% parameters to use in prmflow
Npixel = prmflow.raw.Npixel;
Nslice = prmflow.raw.Nslice;
Nviewprot = prmflow.raw.Nviewprot;
focalspot = prmflow.raw.focalspot;
% gantrytilt = prmflow.raw.gantrytilt;
% gantrytilt = prmflow.protocol.gantrytilt*(pi/180);
focalposition = prmflow.system.focalposition(focalspot, :);
% Nfocal = prmflow.raw.Nfocal;

% detector
detector = prmflow.system.detector;

% fan angles & focal angle(s)
if isfield(prmflow.rebin, 'fanangles')
    fanangles = prmflow.rebin.fanangles;
    focalangle = prmflow.rebin.focalangle;
else
    [fanangles, focalangle] = detpos2fanangles(detector.position, focalposition);
    fanangles = reshape(fanangles, Npixel, Nslice, []);
    prmflow.rebin.fanangles = fanangles;
    prmflow.rebin.focalangle = focalangle;
end

prmflow.rebin = helicalrebinprepare(prmflow.rebin, detector, fanangles, focalangle, Nviewprot);

% viewblock
rebinpipe = prmflow.pipe.(status.nodename);
if isfield(rebinpipe, 'viewblock') && ~isempty(rebinpipe.viewblock)
    viewblock = rebinpipe.viewblock;
else
    viewblock = 48;
end
prmflow.rebin.viewblock = viewblock;
% other rebin prms
prmflow.rebin.Nview = prmflow.raw.Nview;
prmflow.rebin.Nshot = 1;

% low accuracy interpolation
if isfield(prmflow.protocol, 'datablock') && ~isempty(prmflow.protocol.datablock)
    m = nextpow2(prmflow.protocol.datablock);
else
    m = 9;
end
prmflow.rebin.lowaccuracyshift = 23 - m;

% output to recon
prmflow.recon = prmrebin2recon(prmflow.recon, prmflow.rebin);

% pipe line
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpooldata;
    dataflow.buffer.(nodename) = struct();
    dataflow.buffer.(nodename).outpool = struct();
    dataflow.buffer.(nodename).ReadPoint = 1;
    dataflow.buffer.(nodename).WritePoint = 1;
    dataflow.buffer.(nodename).AvailViewindex = -prmflow.rebin.Nextraview(1) + 0;
    dataflow.buffer.(nodename).AvailPoint = dataflow.buffer.(nodename).AvailViewindex;

    dataflow.buffer.(nodename).poolsize = Inf;
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
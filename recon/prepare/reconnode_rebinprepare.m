function [dataflow, prmflow, status] = reconnode_rebinprepare(dataflow, prmflow, status)
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
gantrytilt = prmflow.raw.gantrytilt;
% gantrytilt = prmflow.protocol.gantrytilt*(pi/180);
focalposition = prmflow.system.focalposition(focalspot, :);
Nfocal = prmflow.raw.Nfocal;
scan = prmflow.raw.scan;
angulationcode = prmflow.system.angulationcode;
angulationzero = prmflow.system.angulationzero;

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

if isfield(nodeprm, 'sloperebin')
    prmflow.rebin.issloperebin = nodeprm.sloperebin;
    % suggest to close the sloperebin for 2D-BP and calibration tasks.
else
    % default method
    switch scan
        case 'axial'
            prmflow.rebin.issloperebin = true;
        case 'helical'
            prmflow.rebin.issloperebin = gantrytilt~=0;
            % tilt helical is not a common case, no problem?
        otherwise
            prmflow.rebin.issloperebin = false;
    end
end

% sloperebin-prepare
prmflow.rebin = sloperebinprepare(prmflow.rebin, detector, fanangles, focalangle, Nviewprot, gantrytilt);

% angulation
prmflow.rebin.angulationcode = angulationcode;
prmflow.rebin.angulationzero = angulationzero;
prmflow.rebin.delta_anglecode = angulationcode/Nviewprot;
prmflow.rebin.anglepercode = (pi*2) / angulationcode;
prmflow.rebin.viewangleshift = angulationzero + prmflow.rebin.DFSviewshift;
% I know those are double values

% other rebin prms
prmflow.rebin.Nshot = prmflow.raw.Nshot;
prmflow.rebin.viewnumber = prmflow.raw.viewnumber;
prmflow.rebin.maxviewnumber = prmflow.raw.maxviewnumber;
prmflow.rebin.Nview = prmflow.raw.Nview;
prmflow.rebin.viewread = 0;
prmflow.rebin.viewleft = 0;
% after rebin
prmflow.rebin.NviewOut = floor(prmflow.rebin.Nview / prmflow.rebin.Nfocal);
prmflow.rebin.viewnumberOut = floor(prmflow.rebin.viewnumber / prmflow.rebin.Nfocal);

if isfield(nodeprm, 'lowaccuracyshift')
    % a suggested lowaccuracyshift is 14
    prmflow.rebin = lowaccuracyAzi(prmflow.rebin, nodeprm.lowaccuracyshift);
end

if isfield(nodeprm, 'viewanglealign')
    % on/off flag of viewanglealign
    prmflow.rebin.viewanglealign = nodeprm.viewanglealign;
    % to do this while the frameware can not repeatably keep the accuracy
    % of the gantry rotator's position in shots of scan.
else
    prmflow.rebin.viewanglealign = true;
    % plz make sure if your CT need this patch.
end

% output to recon
prmflow.recon = prmrebin2recon(prmflow.recon, prmflow.rebin);
% othter to output to recon
prmflow.recon.dataclass = prmflow.raw.dataclass;

% % cross shot
% if strcmpi(scan, 'axial') && prmflow.recon.Nshot > 1
%     prmflow.recon.Nshot = prmflow.recon.Nshot + 1;
% end
% the BP nodes will consider that.

% pipe line
if pipeline_onoff
    % pipeline console paramters
    prmflow.pipe.(nodename).pipeline = struct();
    % the rebin correction is H-H.1.G or A.1.G 
    prmflow.pipe.(nodename).pipeline.kernellevel = 1;
    prmflow.pipe.(nodename).pipeline.viewrely = double(prmflow.rebin.Nviewrely);  % int32
    prmflow.pipe.(nodename).pipeline.relystrategy = 'greedy';
    prmflow.pipe.(nodename).pipeline.viewrescale = double([1 Nfocal]);  % int32
    prmflow.pipe.(nodename).pipeline.viewextra = [0 0];
    % ask datasize for next node
    prmflow.pipe.(nodename).pipeline.nextdatasize = double(prmflow.rebin.Nreb * prmflow.rebin.Nslice);
    prmflow.pipe.(nodename).pipeline.nextobjecttype = 'rawdata';
    % ask circulte
    if strcmpi(scan, 'axial')
        prmflow.pipe.(nodename).pipeline.nextcirculte = true;
    end
%     prmflow.pipe.(nodename).pipeline.inputminlimit = min(viewblock, prmflow.raw.Nviewprot);
    
%     % needs to initial the empty rawdata in pipepool
%     nextnode = status.pipeline.(nodename).nextnode;
%     if ~strcmpi(nextnode, 'NULL')
%         dataflow.pipepool.(nextnode) = status.defaultpool;
%         Nps = prmflow.recon.Npixel * prmflow.recon.Nslice;
%         dataflow.pipepool.(nextnode).data.rawdata = zeros(Nps, 0, 'single');
%         % An enmpty rawdata could lay in an error resize step in nodepriostep.m
%     end

    % private buffer
    dataflow.buffer.(nodename) = struct();
    % to save the startangle (of first shot in axial)
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function rebin = lowaccuracyAzi(rebin, lowaccuracyshift)

rebin.lowaccuracyshift = lowaccuracyshift;
rebin.idealfAzi = round(rebin.idealfAzi.*2^lowaccuracyshift).*2^(-lowaccuracyshift);

end


function recon = prmrebin2recon(recon, rebin)
% to output the recon parameters after rebin prepare
%   recon = prmrebin2recon(recon, rebin);
% call it after rebin prepare.

recon.Nshot = rebin.Nshot;
recon.Npixel = rebin.Nreb;
recon.Nslice = rebin.Nslice;
recon.Nfocal = rebin.Nfocal;
recon.delta_d = rebin.delta_d;
recon.delta_z = rebin.delta_z;
recon.delta_view = rebin.delta_view * rebin.Nfocal;
recon.delta_anglecode = rebin.delta_anglecode * rebin.Nfocal;
recon.midchannel = rebin.midU_phi;
recon.SID = rebin.SID;
recon.viewnumber = rebin.viewnumberOut;
recon.Nview = rebin.NviewOut;
recon.maxviewnumber = floor(rebin.maxviewnumber / rebin.Nfocal);
recon.Nviewprot = rebin.Nviewprot / rebin.Nfocal;
recon.gantrytilt = rebin.gantrytilt;
% recon.viewangleshift = rebin.angulationzero + rebin.DFSviewshift; % double

% has not been filtered
recon.filtered = false;

end
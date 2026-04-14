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
FocalMode = prmflow.raw.FocalMode;
gantrytilt = prmflow.raw.gantrytilt;
% gantrytilt = prmflow.protocol.gantrytilt*(pi/180);
focalposition = prmflow.system.focalposition(focalspot, :);
Nfocal = prmflow.raw.Nfocal;
scan = prmflow.raw.scan;
% angulationcode = prmflow.system.angulationcode;
% angulationzero = prmflow.system.angulationzero;
% the prmflow.system shall include the movercode and moverlength to sparse the 'TableEncoder' in rawdata
if isfield(prmflow.system, 'movercode')
    movercode = 2^double(prmflow.system.movercode);
else
    movercode = 2^32;
end
if isfield(prmflow.system, 'moverlength')
    moverlength = double(prmflow.system.moverlength);
else
    moverlength = 2000;
end

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
        case {'helical', 'conveyor'}
            prmflow.rebin.issloperebin = gantrytilt~=0;
            % tilt helical is not a common case, no problem?
        otherwise
            prmflow.rebin.issloperebin = false;
    end
end

% sloperebin-prepare
DFS_onoff = strcmp(FocalMode, 'XDFS');
prmflow.rebin = sloperebinprepare(prmflow.rebin, detector, fanangles, focalangle, Nviewprot, DFS_onoff, gantrytilt);

% angulation
prmflow.rebin.angulationcode = prmflow.raw.angulationcode;
prmflow.rebin.angulationzero = prmflow.raw.angulationzero;
prmflow.rebin.delta_anglecode = double(prmflow.raw.angulationcode)/Nviewprot;
prmflow.rebin.anglepercode = prmflow.raw.anglepercode;
prmflow.rebin.viewangleshift = prmflow.raw.angulationzero + prmflow.rebin.DFSviewshift;
% I know those are double values

% mover
prmflow.rebin.movercode = movercode;
prmflow.rebin.length_movercode = moverlength/movercode;
prmflow.rebin.movershift = prmflow.rebin.DFSviewshift/(pi*2)*Nviewprot * prmflow.rebin.length_movercode;

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

% lowaccuracyshift
if isfield(nodeprm, 'lowaccuracyshift')
    prmflow.rebin.lowaccuracyshift = nodeprm.lowaccuracyshift;
else
    switch scan
        case 'axial'
            prmflow.rebin.lowaccuracyshift = nan;
            % 0 is no low-accuracy
        case {'helical', 'conveyor'}
            prmflow.rebin.lowaccuracyshift = 14;
            % a suggested lowaccuracyshift is 14
        otherwise
            prmflow.rebin.lowaccuracyshift = nan;
    end
end
prmflow.rebin = lowaccuracyAzi(prmflow.rebin);

if isfield(nodeprm, 'viewanglealign')
    % on/off flag of viewanglealign
    prmflow.rebin.viewanglealign = nodeprm.viewanglealign;
    % to do this while the frameware can not repeatably keep the accuracy
    % of the gantry rotator's position in shots of scan.
else
    prmflow.rebin.viewanglealign = true;
    % plz make sure if your CT need this patch.
end

% X-upsampling in rebin
if isfield(nodeprm, 'upsampling') && ~isempty(nodeprm.upsampling)
    prmflow.rebin.upsampling = nodeprm.upsampling;
else
    prmflow.rebin.upsampling = false;
end
if prmflow.rebin.upsampling
    if isfield(nodeprm, 'upsampgamma') && ~isempty(nodeprm.upsampgamma)
        prmflow.rebin.upsampgamma = single(nodeprm.upsampgamma);
    else
        prmflow.rebin.upsampgamma = single([0.7 0.8854]);
    end
else
    prmflow.rebin.upsampgamma = [];
end

% databasis (real or complex)
if ~isfield(prmflow.rebin, 'databasis') || ~isavail(prmflow.rebin.databasis)
    prmflow.rebin.databasis = prmflow.raw.databasis;
end

% output to recon
prmflow.recon = prmrebin2recon(prmflow.recon, prmflow.rebin);
% others to output to recon
prmflow.recon.dataclass = prmflow.raw.dataclass;
prmflow.recon.previewed = false;

if DFS_onoff
    % add 'moverdis' in rawheadfields
    if ~isincell(prmflow.raw.rawheadfields, 'moverdis')
        prmflow.raw.rawheadfields = [prmflow.raw.rawheadfields, {'moverdis'}];
    end
end

% pipe line
if pipeline_onoff
    % pipeline console paramters
    % the rebin correction is H-H.1.G or A.1.G 
    prmflow.pipe.(nodename).pipeline.kernellevel = 1;
    prmflow.pipe.(nodename).pipeline.viewrely = double(prmflow.rebin.Nviewrely);  % int32
    prmflow.pipe.(nodename).pipeline.relystrategy = 'greedy';
    prmflow.pipe.(nodename).pipeline.viewrescale = double([1 Nfocal]);  % int32
    prmflow.pipe.(nodename).pipeline.viewextra = [0 0];
    % ask datasize for next node
    prmflow.pipe.(nodename).pipeline.nextdatasize = double(prmflow.recon.Npixel * prmflow.recon.Nslice);
    prmflow.pipe.(nodename).pipeline.nextobjecttype = 'rawdata';
    % ask circulte
    if strcmpi(scan, 'axial')
        prmflow.pipe.(nodename).pipeline.nextcirculte = true;
    end
    
    % default GPU on
    if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
        prmflow.pipe.(nodename).pipeline.GPUonoff = 1;
    end
    
    % private buffer
    dataflow.buffer.(nodename) = struct();
    % to save the startangle (of first shot in axial)
end

% GPU on/off
prmflow = defaultGPUonoff(prmflow, status, nodename);
% while GPU on
if prmflow.pipe.(nodename).pipeline.GPUonoff > 0
    % put things to GPU
    GPUfields = {'faninterpkern', 'midchannel', 'idealphi', 'Yshift', 'Zgrid', 'idealfAzi', 'lowaccuracyshift'};
    prmflow.rebin = putfieldsinGPU(prmflow.rebin, GPUfields);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function rebin = lowaccuracyAzi(rebin)

if ~isnan(rebin.lowaccuracyshift)
    rebin.idealfAzi = round(rebin.idealfAzi.*2^rebin.lowaccuracyshift).*2^(-rebin.lowaccuracyshift);
end

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
% recon.viewangleshift = rebin.angulationzero + rebin.DFSviewshift; % no longer used
recon.anglepercode = rebin.anglepercode;
recon.viewangleshift = rebin.viewangleshift;
recon.rebinviewrely = ceil(rebin.Nviewrely./rebin.Nfocal);   % it is Rebin.pipeline.viewrely_out
recon.databasis = rebin.databasis;

% has not been filtered
recon.filtered = false;

if rebin.upsampling
    % upsampling
    recon.Npixel = recon.Npixel*2;
    recon.delta_d = recon.delta_d/2;
    recon.midchannel = round(recon.midchannel*4-2)/2;
    recon.upsampled = true;
else
    recon.upsampled = false;
end

end
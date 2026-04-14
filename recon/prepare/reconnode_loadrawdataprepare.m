function [dataflow, prmflow, status] = reconnode_loadrawdataprepare(dataflow, prmflow, status)
% prepare node, read rawdata prepare
% [dataflow, prmflow, status] = reconnode_loadrawdataprepare(dataflow, prmflow, status)

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
pipeline_onoff = status.pipeline.(status.nodename).pipeline_onoff;

% parameters
if isfield(nodeprm, 'maxviewnumber') && ~isempty(nodeprm.maxviewnumber)
    maxviewnumber = nodeprm.maxviewnumber;
else
    maxviewnumber = prmflow.system.maxviewnumber;
end
if isfield(nodeprm, 'restartview') && ~isempty(nodeprm.restartview)
    restartview = nodeprm.restartview;
else
    restartview = prmflow.system.restartview;
end


prmraw = struct();
% angulation
prmraw.angulationcode = prmflow.system.angulationcode;
prmraw.angulationzero = prmflow.system.angulationzero;
prmraw.anglepercode = (pi*2) / double(prmraw.angulationcode);

% .raw from protocol
if isfield(prmflow, 'protocol')
    % Nshot, startshot, endshot
    if isinf(prmflow.protocol.shotnumber) % or shotnumber too large
        % infinite shots
        prmraw.Nshot = ceil(maxviewnumber / prmflow.protocol.viewnumber);
    else
        prmraw.Nshot = prmflow.protocol.shotnumber;
    end
    if isfield(prmflow.protocol, 'startshot')
        prmraw.startshot = prmflow.protocol.startshot;
    else
        prmraw.startshot = 1;
    end
    if isfield(prmflow.protocol, 'endshot')
        prmraw.endshot = prmflow.protocol.endshot;
    else
        prmraw.endshot = prmraw.startshot + prmflow.protocol.shotnumber - 1;
    end
    % recount shotnum
    prmraw.Nshot = min(prmraw.endshot-prmraw.startshot+1, prmraw.Nshot);
    % viewpershot, the protocol.viewnumber is(are) the view-number per shot, but not the total number of the views.
    if length(prmflow.protocol.viewnumber) == 1 
        if ~isinf(prmflow.protocol.shotnumber)
            prmraw.viewpershot = repmat(prmflow.protocol.viewnumber, 1, prmraw.Nshot);
        else
            prmraw.viewpershot = prmflow.protocol.viewnumber;
        end
        startview0 = prmflow.protocol.viewnumber * (prmraw.startshot-1);
    else
        prmraw.viewpershot = prmflow.protocol.viewnumber(prmraw.startshot : end);
        startview0 = sum(prmflow.protocol.viewnumber(1 : prmraw.startshot-1));
    end
    % I know the prmflow.protocol.viewnumber is the view number per shot for axial, and for helical only one shot once.
    % viewnumber, the prmraw.viewnumber is the total number of the views which could not equal to the protocol.viewnumber
    if isinf(prmflow.protocol.shotnumber)
        prmraw.viewnumber = inf;
    else
        prmraw.viewnumber = sum(prmraw.viewpershot);
    end
    % startview, the protocol.startview is the startview to the start shot
    if isfield(prmflow.protocol, 'startview')
        prmraw.startview = prmflow.protocol.startview;
        % fix prmraw.viewnumber
        prmraw.viewnumber = prmraw.viewnumber - prmflow.protocol.startview + 1;
    else
        prmraw.startview = 1;
    end
    % filereadpoint (to be saved in buffer)
    filereadpoint = startview0 + prmraw.startview;
    % datablock onoff
    if pipeline_onoff && ~isempty(prmflow.protocol.datablock)
        prmraw.datablock_onoff = true;
    else
        % do not employ blocked reading in non-pipeline mode
        prmraw.datablock_onoff = false;
    end
    % back up
    prmraw.viewnumber0 = prmraw.viewnumber;
    % Nview is the active viewnumber
    if prmraw.datablock_onoff
        prmraw.Nview = min(prmraw.viewnumber, maxviewnumber);
    else
        prmraw.Nview = prmraw.viewnumber;
    end
    prmraw.maxviewnumber = maxviewnumber;
    prmraw.restartview = restartview;
    % Notes: That,
    % the 'viewnumber0' is the total number of views in a series;
    % the 'viewnumber' is viewnumber0 - viewread that it is the number of views left;
    % in infinite the helical (or mostly converyor), the viewnumber0 and viewnumber are inf;
    % the 'Nview' is the number of active views in this restart period. 

    % Due to the infinite or too large view number is not  accepectable for some recone nodes, we devide them to some restart
    % periods.

    % viewnumber per rotation
    prmraw.Nviewprot = prmflow.protocol.viewperrot;
    % scan
    prmraw.scan = lower(prmflow.protocol.scan);
    % tilt
    prmraw.gantrytilt =  prmflow.protocol.gantrytilt*(pi/180);
    
    % explain focal spot
    focalspot_0x = focalspot20x(prmflow.protocol.focalspot);
    spots = fliplr(dec2bin(focalspot_0x)=='1');
    prmraw.Nfocal = sum(spots);
    prmraw.focalspot = find(spots);
    % NOTE: prmflow.protocol.focalspot is the name of the focalspot mode,
    %       prmflow.raw.focalspot is the index of the focalspot(s).
    prmraw.FocalMode = prmflow.protocol.FocalMode;

    % mover
    prmraw.secperrotation = prmflow.protocol.rotationspeed; % I know the rotationspeed in protocal means second per rotation.
    prmraw.couchspeed = prmflow.protocol.couchspeed;


    % rawdatastyle
    if isfield(prmflow.protocol, 'rawdatastyle') && ~isempty(prmflow.protocol.rawdatastyle)
        prmraw.rawdatastyle = prmflow.protocol.rawdatastyle;
    else
        % default rawdatastyle is '24bit'
        prmraw.rawdatastyle = '24bit';
    end
    % islog2 and dataclass
    switch lower(prmraw.rawdatastyle)
        case '16bit'
            prmraw.islog2 = true;
            prmraw.dataclass = 'real';
        case '24bit'
            prmraw.islog2 = false;
            prmraw.dataclass = 'real';
        case {'float', 'single', 'real'}
            prmraw.islog2 = true;
            prmraw.dataclass = 'real';
        case 'complex'
            prmraw.islog2 = true;
            prmraw.dataclass = 'complex';
        otherwise
            prmraw.islog2 = false;
            prmraw.dataclass = 'real';
    end
    % databasis (only real/complex now)
    if strcmpi(prmraw.dataclass, 'complex')
        prmraw.databasis = 2;
    else
        prmraw.databasis = 1;
    end
    
    % blocked reading
    if prmraw.datablock_onoff
        datablock = prmflow.protocol.datablock;
        Nblock = ceil(prmraw.Nview / datablock);
        prmraw.Nblock = Nblock;
        prmraw.datablocksize = ...
            [repmat(datablock, 1, Nblock-1)  prmraw.Nview - datablock * (Nblock-1)];
    else
        % no blocks
        prmraw.Nblock = 1;
        prmraw.datablocksize = prmraw.Nview;
    end
    % viewread
    prmraw.viewread = 0;
    prmraw.viewleft = 0;
else
    % no protocol data loading (can only load the .mat files)
    prmraw.startshot = 1;
    prmraw.endshot = 1;
    prmraw.Nshot = 1;
    prmraw.viewpershot = nan;
    prmraw.Nblock = 1;
    prmraw.datablocksize = 0;
    % not for recon/cali
end
% ini iblock
prmraw.iblock = 1;
% ini torestart
prmraw.torestart = false;
% rawdata format configure
prmraw.formatcfg = [];
prmraw.formatversion = '';

% copy to prmflow
prmflow.raw = prmraw;

% record Nshot and other in nodeprm
prmflow.pipe.(nodename).Nshot = prmraw.Nshot;
prmflow.pipe.(nodename).viewpershot = prmraw.viewpershot;
prmflow.pipe.(nodename).startshot = prmraw.startshot;
prmflow.pipe.(nodename).endshot = prmraw.endshot;
% prmflow.pipe.(nodename).startview = prmraw.startview;
% I worry some following nodes could change the Nshot

% .raw from calibration tables
if isfield(prmflow.system, 'detector') && ~isempty(prmflow.system.detector)
    prmflow.raw.Nslice = prmflow.system.detector.Nmergedslice;
    prmflow.raw.Npixel = double(prmflow.system.detector.Npixel);
end

% circulate reading
if isfield(nodeprm, 'circulate') && ~isempty(nodeprm.circulate)
    prmflow.raw.circulatereading = nodeprm.circulate;
else
    prmflow.raw.circulatereading = false;
end

% rawhead fields claiming
if ~isfield(prmflow.raw, 'rawheadfields')
    % set default raw head fields
    prmflow.raw.rawheadfields = {'ShotNumber', 'ReadingNumber', 'ShotStart', 'IntegrationTime', 'mA', 'KV', ...
        'AngleEncoder', 'TableEncoder'};
end
if isfield(nodeprm, 'rawheadfields') && ~isempty(nodeprm.rawheadfields)
    newfields = regexp(nodeprm.rawheadfields, '(, +)|(,)', 'split');
    prmflow.raw.rawheadfields = union(prmflow.raw.rawheadfields, newfields);
end
% The 'rawheadfields' has not been employed yet. While it works the subsequent nodes shall can ask new fields in rawhead by
% rewriting the prmflow.raw.rawheadfields like prmflow.raw.rawheadfields = [prmflow.raw.rawheadfields, {newfields}]. The recon
% node readrawdata will carry out those asks.

% private buffer
dataflow.buffer.(nodename) = struct();
% the point in reading file
dataflow.buffer.(nodename).filereadpoint = filereadpoint;
dataflow.buffer.(nodename).ishot = 1;

% initial pipeline
if pipeline_onoff
    % nodename
    nodename = status.nodename;
    % Mostly the nodename is 'loadrawdata' but not 'readrawdata'.
    % node level
    prmflow.pipe.(nodename).pipeline.kernellevel = 1;    % do not copy curr to next
    prmflow.pipe.(nodename).pipeline.inputminlimit = 0;  % need not input data
    % ask objecttype and datasize for next node
    prmflow.pipe.(nodename).pipeline.nextobjecttype = 'rawdata';
    if isfield(prmflow.raw, 'Npixel') && isfield(prmflow.raw, 'Nslice')
        prmflow.pipe.(nodename).pipeline.nextdatasize = prmflow.raw.Npixel*prmflow.raw.Nslice;
    end
    % ask dataclass (not used yet)
    prmflow.pipe.(nodename).pipeline.nextdataclass = prmflow.raw.dataclass;
    % The dataclass should be 'real' or 'complex'. We don't set 'double' or 'single', everything is single yet.

    % input pool
    dataflow.pipepool.(nodename) = status.defaultpool;
    dataflow.pipepool.(nodename).poolsize = double(min(max(prmflow.raw.datablocksize)*4, sum(prmflow.raw.datablocksize)));
    dataflow.pipepool.(nodename).datafields = {'rawdata', 'rawhead'};
    % The 'loadrawdata' is a forced node, the prepare of loadrawdata is also forced. If user configured another node to read
    % rawdata (by calling reconnode_readrawdata), the prepare function reconnode_readrawdataprepare (if exist) will be called to
    % do the prepare for that user configured node. They are different things.
end

% default GPU off
if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
    prmflow.pipe.(nodename).pipeline.GPUonoff = 0;
end

% copy .raw to .prepare
prmflow.prepare = prmflow.raw;
% the prmflow.prepare is specially used by allpurposeprepare and datareloadprepare, not for normal recon jobs.

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end
function [dataflow, prmflow, status] = reconnode_loadrawdataprepare(dataflow, prmflow, status)
% prepare node, read rawdata prepare
% [prmflow, status] = reconnode_loadrawdataprepare(prmflow, status);

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

% pipeline_onoff
pipeline_onoff = status.pipeline.(status.nodename).pipeline_onoff;

prmraw = struct();
% shots
if isfield(prmflow, 'protocol')
    prmraw.Nshot = prmflow.protocol.shotnumber;
    if isfield(prmflow.protocol, 'startshot')
        prmraw.startshot = prmflow.protocol.startshot;
    else
        prmraw.startshot = 1;
    end
    if isfield(prmflow.protocol, 'endshot')
        prmraw.endshot = prmflow.protocol.endshot;
    else
        prmraw.endshot = prmraw.Nshot;
    end
    prmraw.viewpershot = prmflow.protocol.viewnumber;
else
    prmraw.startshot = 1;
    prmraw.endshot = 1;
    prmraw.Nshot = 1;
end
% recount shotnum
prmraw.Nshot = min(prmraw.endshot-prmraw.startshot+1, prmraw.Nshot);
% copy to prmflow
prmflow.raw = prmraw;

% .raw from protocol
if isfield(prmflow, 'protocol')
%     viewnumber = reconcfg.protocol.viewnumber;
    prmflow.raw.Nviewprot = prmflow.protocol.viewperrot;
    % viewnumber
    prmflow.raw.Nview = prmflow.protocol.viewnumber * prmflow.protocol.shotnumber;
    % I know the prmflow.protocol.viewnumber is the view number per shot for axial, and for helical only one shot once.
    % scan
    prmflow.raw.scan = lower(prmflow.protocol.scan);
    % tilt
    prmflow.raw.gantrytilt =  prmflow.protocol.gantrytilt*(pi/180);
    % explain focal spot
    focalspot_0x = focalspot20x(prmflow.protocol.focalspot);
    spots = fliplr(dec2bin(focalspot_0x)=='1');
    prmflow.raw.Nfocal = sum(spots);
    prmflow.raw.focalspot = find(spots);
    % NOTE: prmflow.protocol.focalspot is the name of the focalspot mode,
    %       prmflow.raw.focalspot is the index of the focalspot(s).
    % flag islog2
    if isfield(prmflow.protocol, 'rawdatastyle') && strcmpi(prmflow.protocol.rawdatastyle, '16bit')
        prmflow.raw.islog2 = true;
    else
        prmflow.raw.islog2 = false;
    end
end

% .raw from calibration tables
if isfield(prmflow.system, 'detector') && ~isempty(prmflow.system.detector)
    prmflow.raw.Nslice = prmflow.system.detector.Nmergedslice;
    prmflow.raw.Npixel = double(prmflow.system.detector.Npixel);
end

% initial pipeline
if pipeline_onoff
    % nodename
    nodename = status.nodename;
    % Mostly the nodename is 'loadrawdata' but not 'readrawdata'.

    % private buffer
    if isempty(prmflow.protocol.datablock)
        % no blocks
        dataflow.buffer.(nodename).datablock_onoff = false;
        dataflow.buffer.(nodename).Nblock = 1;
        dataflow.buffer.(nodename).datablocksize = prmflow.protocol.viewnumber * prmflow.protocol.shotnumber;
    else
        dataflow.buffer.(nodename).datablock_onoff = true;
        viewnumber = prmflow.protocol.viewnumber * prmflow.protocol.shotnumber;
        datablock = prmflow.protocol.datablock;
        Nblock = ceil(viewnumber / datablock);
        dataflow.buffer.(nodename).Nblock = Nblock;
        dataflow.buffer.(nodename).datablocksize = ...
            [repmat(datablock, 1, Nblock-1)  viewnumber - datablock * (Nblock-1)];
    end
    dataflow.buffer.(nodename).iblock = 1;
    % to prepare a 'hidden' node of loadrawdata in pipe-line
    dataflow.buffer.(nodename).outputpool = status.defaultpublicpool;
    dataflow.buffer.(nodename).datafields = {'rawdata', 'rawhead'};
    prmflow.pipe.(nodename).level = 1;
    
    % old
    dataflow.buffer.(nodename).pipepool = struct();
    dataflow.buffer.(nodename).pipepool.ReadPoint = 1;
    dataflow.buffer.(nodename).pipepool.WritePoint = 1;
    % hard-codes (we need a readrawdataprepare)
    dataflow.buffer.(nodename).pipepool.poolsize = Inf;
    dataflow.buffer.(nodename).pipepool.warningstage = 2048;
    dataflow.buffer.(nodename).pipepool.recylestrategy = 0;
    % The 'loadrawdata' is a forced node, the prepare of loadrawdata is also forced. If user configured another node to read
    % rawdata (by calling reconnode_readrawdata), the prepare function reconnode_readrawdataprepare (if exist) will be called to
    % do the prepare for that user configured node. They are different things.

    % input pool
    dataflow.pipepool.(nodename) = struct();
    % no input pool fields for the node loadrawdata
end

end
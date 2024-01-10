function [dataflow, prmflow, status] = reconnode_readrawdata(dataflow, prmflow, status)
% recon node to read raw data in recon/cali pipe line
% [dataflow, prmflow, status] = reconnode_readrawdata(dataflow, prmflow, status);
% NOTE: no quick start, if you only want to read a rawdata plz call loadrawdata.m
% this function is a recon pipe line node, but not an I/O function.

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

% not prepared?
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_loadrawdataprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% node name
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

if isempty(prmflow.rawdata)
    % to skip the read rawdata
    [dataflow, prmflow, status] = reconnode_donothing(dataflow, prmflow, status);
    if status.jobdone>0
        status.jobdone = 5;
    end
    return;
elseif ~exist(prmflow.rawdata, 'file')
    status.jobdone = false;
    status.errorcode = 1;
    status.errormsg = '[readrawdata] rawdata file not exist.';
    return;
end

% shots to read
shotnum = prmflow.raw.Nshot;
startshot = prmflow.raw.startshot;
viewpershot = prmflow.raw.viewpershot;
if isfield(dataflow, 'buffer') && dataflow.buffer.(nodename).datablock_onoff
    % to read data by blocks
    datablocksize = dataflow.buffer.(nodename).datablocksize;
    iblock = dataflow.buffer.(nodename).iblock;
    if iblock > length(datablocksize)
        % I think it is a mistake in pipeline
        % pass
        status.jobdone = 3;
        return
    end
    startview = sum(datablocksize(1:iblock-1)) + (startshot-1)*viewpershot + 1;
    viewnum = datablocksize(iblock);
else
    startview = (startshot-1)*viewpershot + 1;
    viewnum = shotnum * viewpershot;
end

% the rawcfg was done?
if isfield(dataflow, 'buffer')
    if isfield(dataflow.buffer.(nodename), 'rawcfg') && ~isempty(dataflow.buffer.(nodename).rawcfg)
        rawcfg = dataflow.buffer.(nodename).rawcfg;
    else
        rawcfg = [];
    end
else
    rawcfg = [];
end

% load raw data
[dataflow, headprm, outcfg] = loadrawdata(dataflow, prmflow.rawdata, prmflow.IOstandard, startview, viewnum, rawcfg);
% dataflow = structmerge(loadrawdata(prmflow.rawdata, prmflow.IOstandard, startview, viewnum, rawcfg), dataflow, 0, 0);

% load offset
if isfield(prmflow, 'offset') && ~isempty(prmflow.offset) && ~isfield(dataflow, 'offset')
    dataflow.offset = loadrawdata(struct(), prmflow.offset, prmflow.IOstandard, 1, inf, rawcfg);
end

% a hard code to fill up viewangle by Angle_encoder
if isfield(dataflow.rawhead, 'Angle_encoder') && isfield(prmflow, 'system')
    dataflow.rawhead.viewangle = (single(dataflow.rawhead.Angle_encoder) - prmflow.system.angulationzero) ...
        ./prmflow.system.angulationcode.*(pi*2);
end
% to fill up the Reading_Number if not exist
if ~isfield(dataflow.rawhead, 'Reading_Number')
    dataflow.rawhead.Reading_Number = startview : startview+viewnum-1;
end
% to fill up the Shot_Number if not exist
if ~isfield(dataflow.rawhead, 'Shot_Number')
    dataflow.rawhead.Shot_Number = ceil(dataflow.rawhead.Reading_Number / viewpershot);
end

% prmflow.rawhead
if ~isfield(prmflow, 'rawhead') || isempty(prmflow.rawhead)
    prmflow.rawhead = headprm;
end

% % copy to recon parameters
% % shot number
% prmflow.recon.Nshot = shotnum;
% % from protocol
% if isfield(prmflow, 'protocol')
% %     viewnumber = reconcfg.protocol.viewnumber;
%     prmflow.recon.Nviewprot = prmflow.protocol.viewperrot;
%     % viewnumber
%     prmflow.recon.Nview = prmflow.protocol.viewnumber * shotnum;
%     % I know the prmflow.protocol.viewnumber is the view number per shot for axial, and for helical only one shot once.
%     % scan
%     prmflow.recon.scan = lower(prmflow.protocol.scan);
%     % tilt
%     prmflow.recon.gantrytilt =  prmflow.protocol.gantrytilt*(pi/180);
%     % explain focal spot
%     focalspot_0x = focalspot20x(prmflow.protocol.focalspot);
%     spots = fliplr(dec2bin(focalspot_0x)=='1');
%     prmflow.recon.Nfocal = sum(spots);
%     prmflow.recon.focalspot = find(spots);
%     % NOTE: prmflow.protocol.focalspot is the name of the focalspot mode,
%     %       prmflow.recon.focalspot is the index of the focalspot(s).
%     
% end

% datablock
if isfield(dataflow, 'buffer') && dataflow.buffer.(nodename).datablock_onoff
    % contol the looping of data blocks
    iblock = dataflow.buffer.(nodename).iblock;
    Nblock = dataflow.buffer.(nodename).Nblock;
    if iblock == 1
        dataflow.buffer.(nodename).cfg = outcfg;
    end
    if iblock < Nblock
        status.jobdone = 2;
    else
%         status.pipeline.(nodename).sleeping = true;
        status.jobdone = 1;
    end
    % iblock++
    dataflow.buffer.(nodename).iblock = iblock + 1;
    % WritePoint for the 'hidden' node
    if isfield(dataflow, 'buffer')
        dataflow.buffer.(nodename).pipepool.WritePoint = dataflow.buffer.(nodename).pipepool.WritePoint + viewnum;
    end
else
    status.jobdone = true;
end
% The 'loadrawdata' is a special node in pipe-line. It reads data from file to dataflow.rawdata (and dataflow.head), but they
% are not in the pipeline pool for that we need a 'hidden' node to copy the .rawdata and .head to the next node's input pool.
% The status.pipepool.loadrawdata.WritePoint and ReadPoint are pointing to the dataflow.rawdata but not the loadrawdata's
% input pool or the file handle. 
% The position of the file handle to start reading is saved in dataflow.buffer.loadrawdata ('iblock' somehow); 
% The dataflow.pipepool.loadrawdata is normally empty due to it is the first node, no nodes send data to it;
% That 'hidden' node will 'read' data from dataflow.rawdata and copy to next node's input pool.
if pipeline_onoff && ~dataflow.buffer.(nodename).datablock_onoff
    % What? Shouldn't we forbid this behavior?
    dataflow.buffer.(nodename).pipepool.WritePoint = dataflow.buffer.(nodename).pipepool.WritePoint + viewnum;
end

% pipeline hidden node
if pipeline_onoff
    nextnode = status.pipeline.(nodename).nextnode;
    % the current pool of the 'hidden' node is not in status.pipepool, but in private buffer
    currpool = dataflow.buffer.(nodename).pipepool;
    if ~isempty(nextnode)
        statusnext = status.pipepool.(nextnode);
        % data size to be written in the next pool
        writenum = min(viewnum, min(statusnext.WriteEnd, statusnext.poolsize) - statusnext.WritePoint + 1);
        % copy rawdata and rawhead to next pool
        1;
        [dataflow.pipepool.(nextnode), writenum] = pooldatacopy(dataflow, dataflow.pipepool.(nextnode), ...
            currpool.ReadPoint, statusnext.WritePoint, writenum, {'rawdata', 'rawhead'}, true);

        % move next pool's write point
        status.pipepool.(nextnode).WritePoint = status.pipepool.(nextnode).WritePoint + writenum;
        % Maybe we shall let pipeline consol doing that.
        % move current pool's read point
        dataflow.buffer.(nodename).pipepool.ReadPoint = dataflow.buffer.(nodename).pipepool.ReadPoint + writenum;
        % recycle
        currpool = dataflow.buffer.(nodename).pipepool;
        % datalength
        datalength = currpool.WritePoint - currpool.ReadPoint;
        switch currpool.recylestrategy
            case 0
                recycle_onoff = false;
            case 1
                % recycle
                recycle_onoff = true;
            case 2
                % recycle?
                recycle_onoff = currpool.WritePoint > currpool.warningstage;
            otherwise
                status.jobdone = false;
                status.errorcode = 901;
                status.errormsg = sprintf('readrawdata''s hidden node error, unkown recylestrategy %d!', currpool.recylestrategy);
                return
        end
        if recycle_onoff && datalength>0
            dataflow.rawdata(:, 1:datalength) = ...
                dataflow.rawdata(:, currpool.ReadPoint : currpool.WritePoint-1);
            headfields = fieldnames(dataflow.rawhead);
            for ii = 1:length(headfields)
                dataflow.rawhead.(headfields{ii})(:, 1:datalength) = ...
                    dataflow.rawhead.(headfields{ii})(:, currpool.ReadPoint : currpool.WritePoint-1);
            end
            dataflow.buffer.(nodename).pipepool.ReadPoint = 1;
            dataflow.buffer.(nodename).pipepool.WritePoint = datalength + 1;
        end
        
    end 
end
% done

status.errorcode = 0;
status.errormsg = [];

end


function [dataflow, headprm, outcfg] = loadrawdata(dataflow, filename, IOstandard, startview, viewnum, rawcfg)
% load rawdata (or offset) from the file filename

% protocol = [];
outcfg = [];
headprm = [];
[~, ~, fileEXT] = fileparts(filename);
switch lower(fileEXT)
    case {'.raw', '.bin'}
        % tmp code
        % raw = loaddata(filename, IOstandard);
        % .raw should have a single code to read
        % yep, let'd do it
        if isempty(rawcfg)
            rawcfg = readcfgfile(cfgmatchrule(filename, IOstandard));
        end
        rawcfg.number = viewnum;
        fid = fopen(filename, 'r');
        [raw, outcfg] = sparsepack(fid, rawcfg, startview-1);
        fclose(fid);
        % data flow
        [rawhead, rawdata, headprm] = raw2dataflow(raw);
        dataflow = rawdatamerge(dataflow, rawhead, rawdata);
    case '.mat'
        % load mat
        raw = load(filename);
        % data flow
        if isfield(raw, 'rawhead') && isfield(raw, 'rawdata')
            dataflow = rawdatamerge(dataflow, raw.rawhead, raw.rawdata);
            % but we can not select the views after it has been merged in dataflow
            if isfield(raw, 'headprm')
                headprm = raw.headprm;
            end
        else
            tmpfield = fieldnames(raw);
            [rawhead, rawdata, headprm] = raw2dataflow(raw.(tmpfield{1}), startview, viewnum);
            dataflow = rawdatamerge(dataflow, rawhead, rawdata);
        end
    case '.pd'
        % external IO of .pd
        [raw, protocol] = CRIS2dataflow(filename, startview, viewnum);
        dataflow = rawdatamerge(dataflow, raw.rawhead, raw.rawdata, raw.offset);
        % We need a function to tranlate the CRIS protocol to headprm
        headprm = protocol;
    otherwise
        if  regexp(fileEXT, '[.]m\d+$')
            % the file like .m01, .m02
            raw = load(filename, '-mat');
            tmpfield = fieldnames(raw);
            [rawhead, rawdata, headprm] = raw2dataflow(raw.(tmpfield{1}));
            % In this case the datablock is fixed (by file size)
            dataflow = rawdatamerge(dataflow, rawhead, rawdata);
        else
            error('Unknown rawdata ext: %s', fileEXT);
        end
end

end


function [rawhead, rawdata, headprm] = raw2dataflow(raw, startview, viewnum)
% raw to dataflow

% current shot(s)
if nargin > 2
    Nraw = size(raw(:),1);
    endview = min(startview + viewnum - 1, Nraw);
    raw = raw(startview : endview);
elseif nargin > 1
    raw = raw(startview : end);
end

rawhead.Angle_encoder = [raw.Angle_encoder];
rawhead.Reading_Number = [raw.Reading_Number];
rawhead.Integration_Time = [raw.Integration_Time];
rawhead.Shot_Number = [raw.Shot_Number];
% rawhead.Time_Stamp = [raw.Time_Stamp];
rawhead.mA = single([raw.mA]);
rawhead.KV = single([raw.KV]);
rawdata = single([raw.Raw_Data]);

headprm.Package_Version = raw(1).Package_Version;
headprm.Series_Number = raw(1).Series_Number;
headprm.Start_Slice = raw(1).Start_Slice;
headprm.End_Slice = raw(1).End_Slice;
headprm.Slice_mergescale = raw(1).Slice_mergescale;
headprm.Slice_Number = raw(1).Slice_Number;
headprm.Raw_Data_Size = raw(1).Raw_Data_Size;

end


function dataflow = rawdatamerge(dataflow, rawhead, rawdata, offset)

% rawhead
if isfield(dataflow, 'rawhead')
    headfields = fieldnames(rawhead);
    Nhf = size(headfields(:),1);
    for ihf = 1:Nhf
        hfield_ii = headfields{ihf};
        if isfield(dataflow.rawhead, hfield_ii)
            dataflow.rawhead.(hfield_ii) = [dataflow.rawhead.(hfield_ii) rawhead.(hfield_ii)(:)'];
        else
            dataflow.rawhead.(hfield_ii) = rawhead.(hfield_ii)(:)';
        end
    end
else
    dataflow.rawhead = rawhead;
end

% rawdata
if isfield(dataflow, 'rawdata')
    dataflow.rawdata = [dataflow.rawdata rawdata];
else
    dataflow.rawdata = rawdata;
end

% offset
if nargin > 3 && isavail(offset)
    if isfield(dataflow, 'offset')
        rawdatamerge(dataflow.offset, offset.rawhead, offset.rawdata)
    else
        dataflow.offset = offset;
    end
end

end


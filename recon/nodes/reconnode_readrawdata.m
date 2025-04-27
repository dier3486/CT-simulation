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

if isempty(prmflow.protocol.rawdata)
    % to skip the read rawdata
    status.jobdone = 5;
    return;
elseif ~exist(prmflow.protocol.rawdata, 'file') && ~strcmpi(prmflow.protocol.rawdata, 'fake')
    status.jobdone = false;
    status.errorcode = 1;
    status.errormsg = '[readrawdata] rawdata file not exist.';
    return;
end

% shots to read
shotnum = prmflow.raw.Nshot;
startshot = prmflow.raw.startshot;
viewpershot = prmflow.raw.viewpershot;
% if prmflow.raw.datablock_onoff was prepared

% to read data by blocks (or one block)
datablocksize = prmflow.raw.datablocksize;
iblock = prmflow.raw.iblock;
Nblock = prmflow.raw.Nblock;
if iblock <= Nblock
    readingdata = true;
    % startview, viewnum
%     if length(viewpershot) == 1
%         startview = sum(datablocksize(1:iblock-1)) + (startshot-1)*viewpershot + 1;
%     else
%         startview = sum(datablocksize(1:iblock-1)) + sum(viewpershot(1 : startshot-1)) + 1;
%     end
    viewnum = datablocksize(iblock);
    startview = dataflow.buffer.(nodename).filereadpoint;
    % record read view number
    prmflow.raw.viewread = prmflow.raw.viewread + viewnum;
    % iblock++
    prmflow.raw.iblock = iblock + 1;
else
    % to pass
    readingdata = false;
end

% the rawdata formatcfg was done?
if isfield(prmflow.raw, 'formatcfg')
    rawcfg = prmflow.raw.formatcfg;
else
    rawcfg = [];
end
% the formatcfg is the file format configure struct of in the reading raw data.

% to read the file in circulate
circulatereading = prmflow.raw.circulatereading;

if readingdata
    % load raw data
    if pipeline_onoff
        [dataflow.pipepool.(nodename), dataflow.pipepool.(nodename).data, offset, headprm, outcfg] = ...
            loadrawdata(dataflow.pipepool.(nodename), dataflow.pipepool.(nodename).data, prmflow.protocol.rawdata, ...
            prmflow, startview, viewnum, rawcfg, circulatereading);
    else
        [~, dataflow, offset, headprm, outcfg] = loadrawdata([], dataflow, prmflow.protocol.rawdata, ...
            prmflow, startview, viewnum, rawcfg, circulatereading);
    end

    % load offset
    if ~isfield(dataflow, 'offset')
        if isfield(prmflow, 'offset') && ~isempty(prmflow.offset)
            % read offset from prmflow.offset
            [~, dataflow.offset] = loadrawdata([], struct(), prmflow.offset, ...
                prmflow, 1, inf, rawcfg);
        else
            % it was read from rawdata file or it is empty
            dataflow.offset = offset;
        end
    end
    
    % file read point
    dataflow.buffer.(nodename).filereadpoint = dataflow.buffer.(nodename).filereadpoint + viewnum;
    if circulatereading
        dataflow.buffer.(nodename).filereadpoint = mod(dataflow.buffer.(nodename).filereadpoint, outcfg.filelength);
    end
end

% record rawhead
if ~isfield(prmflow.raw, 'rawhead') || isempty(prmflow.raw.rawhead)
    prmflow.raw.rawhead = headprm;
end
% record outcfg
if ~isfield(prmflow.raw, 'formatcfg')
    prmflow.raw.formatcfg = outcfg;
end

% to restart
if prmflow.raw.iblock >= Nblock
    % I know the prmflow.raw.iblock was +1
    viewread = double(prmflow.raw.viewread);
    if prmflow.raw.iblock == Nblock
        viewread = viewread + prmflow.raw.datablocksize(Nblock);
    end
    status.torestart = viewread < prmflow.raw.viewnumber;
else
    status.torestart = false;
end

% contol the looping of data blocks
if prmflow.raw.iblock <= prmflow.raw.Nblock || status.torestart
    % keep waking
    status.jobdone = 2;
else
    status.jobdone = 1;
end
% I know while the prmflow.raw.datablock_onoff is false, the Nblock=1.

% The 'loadrawdata' is a special node in pipe-line. It reads data from file to dataflow.rawdata (and dataflow.head), but they
% are not in the pipeline pool for that we need a 'hidden' node to copy the .rawdata and .head to the next node's input pool.
% The dataflow.pipepool.loadrawdata.WritePoint and ReadPoint are pointing to the dataflow.rawdata but not the loadrawdata's
% input pool or the file handle. 
% The position of the file handle to start reading is saved in dataflow.buffer.loadrawdata ('iblock' somehow); 
% The dataflow.pipepool.loadrawdata is normally empty due to it is the first node, no nodes send data to it;
% That 'hidden' node will 'read' data from dataflow.rawdata and copy to next node's input pool.

% pipeline hidden node
if pipeline_onoff
    % data can output
    AvailNumber = dataflow.pipepool.(nodename).AvailPoint - dataflow.pipepool.(nodename).ReadPoint + 1;
    % nextnode
    nextnode = status.pipeline.(nodename).nextnode;
    if ~isempty(dataflow.pipepool.(nextnode)) && ~dataflow.pipepool.(nextnode).WriteStuck
        % check if shot start
        isshotstart = dataflow.pipepool.(nodename).ReadPoint == dataflow.pipepool.(nodename).ReadStart || ...
            dataflow.pipepool.(nodename).ReadPoint > dataflow.pipepool.(nodename).ReadEnd;
%         if ~isavail(dataflow.pipepool.(nextnode).WriteEnd) && ~dataflow.pipepool.(nextnode).circulatemode
        % set WriteEnd for next pool
        if isshotstart
            if ~isfield(dataflow.buffer.(nodename), 'ishot')
                dataflow.buffer.(nodename).ishot = 0;
            end
            % view number of current shot (Nviewcurrshot)
            ishot = dataflow.buffer.(nodename).ishot + startshot;
            if length(viewpershot) == 1
                Nviewcurrshot = double(viewpershot);
            else
                Nviewcurrshot = double(viewpershot(ishot));
            end
            % reset currpool's ReadStart and ReadEnd
            dataflow.pipepool.(nodename).ReadStart = dataflow.pipepool.(nodename).ReadPoint;
            dataflow.pipepool.(nodename).ReadEnd = dataflow.pipepool.(nodename).ReadStart + Nviewcurrshot - 1;
            % reset nextpool's WriteStart and WriteEnd
            dataflow.pipepool.(nextnode).WriteStart = dataflow.pipepool.(nextnode).WritePoint;
            dataflow.pipepool.(nextnode).WriteEnd = dataflow.pipepool.(nextnode).WriteStart + Nviewcurrshot - 1;
            % reset nextpool's ReadPoint
            dataflow.pipepool.(nextnode).ReadPoint = dataflow.pipepool.(nextnode).WritePoint;
            dataflow.pipepool.(nextnode).ReadStart = dataflow.pipepool.(nextnode).ReadPoint;
            dataflow.pipepool.(nextnode).ReadEnd = dataflow.pipepool.(nextnode).ReadStart + Nviewcurrshot - 1;
            % reset nextpool's AvialPoint
            dataflow.pipepool.(nextnode).AvailPoint = dataflow.pipepool.(nextnode).ReadPoint-1;
            % close the shotstart
            dataflow.pipepool.(nextnode).isshotstart = false;
            % check poolsize
            if dataflow.pipepool.(nextnode).circulatemode
                if dataflow.pipepool.(nextnode).poolsize ~= Nviewcurrshot
                    % should be some buffer re-malloc here.
                    dataflow.pipepool.(nextnode).poolsize = Nviewcurrshot; 
                end
            end
            % Here is a check point for C++ codes to new or re-malloc the buffers.
            
            % ishot++ 
            dataflow.buffer.(nodename).ishot = dataflow.buffer.(nodename).ishot + 1;
        end
        % nextpoolleft
        nextpoolleft = poolpspaceleft(dataflow.pipepool.(nextnode));
        % I know, after the WritePoint of the next pool reaching the WriteEnd whose WriteStuck will be locked (by
        % movepointsaftercopy.m) and will only be unlocked by poolrecycle.m after the AvailNumber consumed to 0.
        
        % carrynode
        if  dataflow.pipepool.(nextnode)(1).iscarried
            carrynode = dataflow.pipepool.(nextnode)(1).carrynode;
        else
            carrynode = nextnode;
        end

        % copy rawdata to next pool
        writenum = min(nextpoolleft, AvailNumber);
        dataflow.pipepool.(carrynode)(1).data = pooldatacopy(dataflow.pipepool.(nodename), dataflow.pipepool.(nodename).data, ...
            dataflow.pipepool.(nextnode)(1), dataflow.pipepool.(carrynode)(1).data, writenum, [], true);
        % We shall use the carrynode to replace the nextnode while calling the .data.
    elseif isempty(dataflow.pipepool.(nextnode))
        % next node is NULL
        % to return the writenum
        writenum = AvailNumber;
    else
        % nextpool WriteStucked
        writenum = 0;
    end
    % to return the readnumber/writenumber
    status.currentjob.pipeline.readnumber = writenum;
    status.currentjob.pipeline.writenumber = writenum;
    status.currentjob.pipeline.newAvail = writenum;
    status.currentjob.pipeline.Nexpand = 0;
    
    if writenum == 0 && AvailNumber>0
        status.jobdone = 6;
    elseif writenum < AvailNumber
        status.jobdone = 2;
    end
    % I know, else while writenum=AvailNumber>0, status.jobdone=1 if iblock > Nblock, or status.jobdone=2

    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);

end
% done

end


function [currpool, currdata, offset, headprm, outcfg] = loadrawdata(currpool, currdata, filename, ...
    prmflow, startview, viewnum, rawcfg, circulatereading)
% load rawdata (or offset) from the file filename to dataflow (or currdata)
% Note that not all the tags in the rawdata can be read to the dataflow, see function raw2dataflow, only a limited part of them
% will be extracted to support most cali and recon tasks.

if nargin < 8
    circulatereading = false;
end

if isempty(currpool)
    writePoint = [];
else
    writePoint = currpool.WritePoint;
end

offset = [];
outcfg = [];
headprm = [];
IOstandard = prmflow.IOstandard;

if strcmpi(filename, 'fake')
    fileEXT = 'fake';
else
    [~, ~, fileEXT] = fileparts(filename);
end

switch lower(fileEXT)
    case {'.raw', '.bin'}
        % tmp code
        % raw = loaddata(filename, IOstandard);
        % .raw should have a single code to read
        % yep, let'd do it
        if isempty(rawcfg)
            rawcfg = readcfgfile(cfgmatchrule(filename, IOstandard));
        end
        readingnumber = viewnum;
        rawcfg.number = readingnumber;
        viewskip = startview-1;
        if circulatereading && isfield(rawcfg, 'filelength')
            viewskip = mod(viewskip, rawcfg.filelength);
        end
        fid = fopen(filename, 'r');
        [raw, outcfg] = sparsepack(fid, rawcfg, viewskip);
        viewnum = outcfg.number;
        while circulatereading && outcfg.number < readingnumber
            readingnumber = readingnumber - outcfg.number;
            outcfg.number = readingnumber;
            fseek(fid, 0, 'bof');
            [raw1, outcfg] = sparsepack(fid, outcfg, 0);
            raw = cat(2, raw, raw1);
            viewnum = viewnum + outcfg.number;
        end
        fclose(fid);
        % data flow
        [rawhead, rawdata, headprm] = raw2dataflow(raw);
        1;
    case '.mat'
        % load mat
        raw = load(filename);
        % data flow
        if isfield(raw, 'rawhead') && isfield(raw, 'rawdata')
            rawhead = raw.rawhead;
            rawdata = raw.rawdata;
            % but we can not select the views after it has been merged in dataflow
            if isfield(raw, 'headprm')
                headprm = raw.headprm;
            end
            if isfield(raw, 'offset')
                offset = raw.offset;
            end
        else
            tmpfield = fieldnames(raw);
            [rawhead, rawdata, headprm] = raw2dataflow(raw.(tmpfield{1}), startview, viewnum);
        end
    case '.pd'
        % external IO of .pd
        [rawdata, rawhead, offset, headprm] = CRIS2dataflow(filename, startview, viewnum);
        % We may need a function to tranlate the CRIS protocol to headprm
    case 'fake'
        % fake rawdata
        rawdata = ones(1, viewnum, 'single');
        rawhead = struct();
        rawhead.OrigReadingNumber = single(startview : startview + viewnum - 1);
    otherwise
        if  regexp(fileEXT, '[.]m\d+$')
            % the file like .m01, .m02
            raw = load(filename, '-mat');
            tmpfield = fieldnames(raw);
            [rawhead, rawdata, headprm] = raw2dataflow(raw.(tmpfield{1}));
            % In this case the datablock is fixed (by file size)
        else
            error('Unknown rawdata ext: %s', fileEXT);
        end
end

% fill up default rawhead
rawhead = defaultrawhead(rawhead, prmflow.system, startview, viewnum, prmflow.raw.viewpershot);

% write rawdata to currdata
currdata = rawdatamerge(currdata, rawhead, rawdata, writePoint);

% move the points in currpool
if ~isempty(currpool)
    currpool.WritePoint = currpool.WritePoint + double(viewnum);
    currpool.AvailPoint = currpool.AvailPoint + double(viewnum);
end

end

function [rawhead, rawdata, headprm] = raw2dataflow(raw, startview, viewnum)
% raw to dataflow structure format

% raw shall be in size 1*Nview
raw = raw(:)';

% current shot(s)
if nargin > 2
    Nview = size(raw, 2);
    endview = min(startview + viewnum - 1, Nview);
    raw = raw(startview : endview);
elseif nargin > 1
    raw = raw(startview : end);
end

% rawhead
rawhead.Angle_encoder = [raw.Angle_encoder];
rawhead.Reading_Number = [raw.Reading_Number];
rawhead.Integration_Time = [raw.Integration_Time];
rawhead.Shot_Number = [raw.Shot_Number];
% rawhead.Time_Stamp = [raw.Time_Stamp];
rawhead.mA = single([raw.mA]);
rawhead.KV = single([raw.KV]);
rawhead.Table_encoder = [raw.Table_encoder];
if isfield(raw, 'Table_gear')
    rawhead.Table_gear = [raw.Table_gear];
end
% We may have a list to configure in reading those tags to rawhead, which could be depending on protocol.

% rawdata
rawdata = single([raw.Raw_Data]);
% I know all the fields of rawhead and rawdata are in size n*Nview

if ~isempty(raw)
    headprm.Package_Version = raw(1).Package_Version;
    headprm.Series_Number = raw(1).Series_Number;
    headprm.Start_Slice = raw(1).Start_Slice;
    headprm.End_Slice = raw(1).End_Slice;
    headprm.Slice_mergescale = raw(1).Slice_mergescale;
    headprm.Slice_Number = raw(1).Slice_Number;
    headprm.Raw_Data_Size = raw(1).Raw_Data_Size;
else
    headprm = [];
end

end

function dataflow = rawdatamerge(dataflow, rawhead, rawdata, writepoint)
% to write the rawdata to dataflow (or pool) 
if nargin < 4
    writepoint = [];
end

% rawhead
if isfield(dataflow, 'rawhead')
    headfields = fieldnames(rawhead);
    Nhf = size(headfields(:),1);
    for ihf = 1:Nhf
        hfield_ii = headfields{ihf};
        if isfield(dataflow.rawhead, hfield_ii)
            if isempty(writepoint)
                dataflow.rawhead.(hfield_ii) = [dataflow.rawhead.(hfield_ii) rawhead.(hfield_ii)];
            else
                n = size(rawhead.(hfield_ii), 2);
                dataflow.rawhead.(hfield_ii)(:, writepoint : writepoint+n-1) = rawhead.(hfield_ii);
            end
        else
            dataflow.rawhead.(hfield_ii) = rawhead.(hfield_ii);
        end
    end
else
    dataflow.rawhead = rawhead;
end

% rawdata
if isfield(dataflow, 'rawdata')
    if isempty(writepoint)
        dataflow.rawdata = [dataflow.rawdata rawdata];
    else
        n = size(rawdata, 2);
        dataflow.rawdata(:, writepoint : writepoint+n-1) = rawdata;
    end
else
    dataflow.rawdata = rawdata;
end

% % offset
% if nargin > 4 && isavail(offset)
%     if isfield(dataflow, 'offset')
%         rawdatamerge(dataflow.offset, offset.rawhead, offset.rawdata);
%     else
%         dataflow.offset = offset;
%     end
% end

end

function rawhead = defaultrawhead(rawhead, system, startview, viewnum, viewpershot)
% hard coded default rawhead
% to fill up the missing fields and align the data class

% to fill up viewangle by Angle_encoder
if isfield(rawhead, 'Angle_encoder')
    rawhead.Angle_encoder = uint32(rawhead.Angle_encoder);
    if ~isfield(rawhead, 'viewangle') && isfield(system, 'angulationcode')
        rawhead.viewangle = (single(rawhead.Angle_encoder) - system.angulationzero) ...
            ./system.angulationcode.*(pi*2);
    end
end

% viewindex and shotindex
viewindex = double(startview) : double(startview+viewnum)-1;
viewpershot_cum = [0; cumsum(viewpershot(:))];
shotindex = sum(viewindex > viewpershot_cum, 1);

% to fill up the Reading_Number if not exist
if ~isfield(rawhead, 'Reading_Number')
    rawhead.Reading_Number = uint32(viewindex - viewpershot_cum(shotindex)');
else
    rawhead.Reading_Number = uint32(rawhead.Reading_Number);
end
% to fill up the Shot_Number if not exist
if ~isfield(rawhead, 'Shot_Number')
    rawhead.Shot_Number = uint16(shotindex);
else
    rawhead.Shot_Number = uint16(rawhead.Shot_Number);
end
% to fill up the Shot_Start if not exist
if ~isfield(rawhead, 'Shot_Start')
    currShotStart = viewindex - viewpershot_cum;
    rawhead.Shot_Start = int16(any(currShotStart==1, 1)) - int16(any(currShotStart==0, 1));
else
    rawhead.Shot_Start = int16(rawhead.Shot_Start);
end
% the Shot_Start is demanded to be 1 for the first view of a shot, -1 for the end view of a shot and 0 for the others. 
% Note: in pipeline console we will not use the rawhead.Shot_Start to judge the shot start/end. Actually we never use any
% information from the the dataflow to govern the pipeline accesses.
end

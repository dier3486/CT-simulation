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

% pipelineOnoff
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

% to restart the pipeline
if prmflow.raw.torestart
    status.torestart = true;
    prmflow.raw.torestart = false;
    status.jobdone = 7;     % keep waking
    return;
end

% shots to read
% viewpershot = prmflow.raw.viewpershot;
viewpershot = prmflow.pipe.(nodename).viewpershot;
startview = prmflow.raw.startview;

% to read data by blocks (or one block)
datablocksize = prmflow.raw.datablocksize;
iblock = prmflow.raw.iblock;
Nblock = prmflow.raw.Nblock;
viewrestart = prmflow.raw.maxviewnumber - prmflow.raw.restartview;

if iblock <= Nblock
    readingdata = true;
else
    % to pass
    readingdata = false;
end

if readingdata
    % view to read
    viewnum = datablocksize(iblock);
    % view read
    viewread = prmflow.raw.viewread;
    if iblock > 1 && viewnum + viewread > viewrestart
        % to restart
        status.torestart = true;
        status.jobdone = 7;     % keep waking
        return;
    end
    filereadpoint = dataflow.buffer.(nodename).filereadpoint;
    % record read view number
    prmflow.raw.viewread = prmflow.raw.viewread + viewnum;
    % iblock++
    prmflow.raw.iblock = iblock + 1;
end

% % the rawdata formatcfg was done?
% if isfield(prmflow.raw, 'formatcfg')
%     rawcfg = prmflow.raw.formatcfg;
% else
%     rawcfg = [];
% end
% the formatcfg is the file format configure struct of in the reading raw data.

% to read the file in circulate
circulatereading = prmflow.raw.circulatereading;

if readingdata
    % load raw data
    if pipeline_onoff
        [dataflow.pipepool.(nodename), dataflow.pipepool.(nodename).data, offset, headprm, outcfg] = ...
            loadrawdata(dataflow.pipepool.(nodename), dataflow.pipepool.(nodename).data, prmflow.protocol.rawdata, ...
            prmflow, filereadpoint, viewnum, circulatereading);
    else
        [~, dataflow, offset, headprm, outcfg] = loadrawdata([], dataflow, prmflow.protocol.rawdata, ...
            prmflow, filereadpoint, viewnum, circulatereading);
    end
    % TBC: The loadrawdata shall check the prmflow.raw.rawheadfields to fill up the missing rawhead fields in rawdata.

    % load offset
    if ~isfield(dataflow, 'offset')
        if isfield(prmflow, 'offset') && ~isempty(prmflow.offset)
            % read offset from prmflow.offset
            [~, dataflow.offset] = loadrawdata([], struct(), prmflow.offset, ...
                prmflow, 1, inf);
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
if ~isfield(prmflow.raw, 'headprm') || isempty(prmflow.raw.headprm)
    prmflow.raw.headprm = headprm;
end
% record outcfg
if ~isfield(prmflow.raw, 'formatcfg')
    prmflow.raw.formatcfg = outcfg;
    if isfield(headprm, 'versionflag')
        prmflow.raw.formatversion = headprm.versionflag;
    end
end

% contol the looping of data blocks
if prmflow.raw.iblock <= prmflow.raw.Nblock
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
    if ~isempty(dataflow.pipepool.(nextnode)) && ~dataflow.pipepool.(nextnode)(1).WriteStuck
        % check if shot start
        isshotstart = dataflow.pipepool.(nodename).ReadPoint == dataflow.pipepool.(nodename).ReadStart || ...
            dataflow.pipepool.(nodename).ReadPoint > dataflow.pipepool.(nodename).ReadEnd;
        status.currentjob.pipeline.isshotstart = isshotstart && dataflow.pipepool.(nextnode)(1).isshotstart;
        % set WriteEnd for next pool
        if isshotstart
            ishot = dataflow.buffer.(nodename).ishot;
            if length(viewpershot) == 1
                Nviewcurrshot = double(viewpershot);
            else
                Nviewcurrshot = double(viewpershot(ishot));
            end
            if ishot == 1
                Nviewcurrshot = Nviewcurrshot - startview + 1;
            end
            % reset currpool's ReadStart and ReadEnd
            dataflow.pipepool.(nodename).ReadStart = dataflow.pipepool.(nodename).ReadPoint;
            dataflow.pipepool.(nodename).ReadEnd = dataflow.pipepool.(nodename).ReadStart + Nviewcurrshot - 1;
            % reset nextpool's WriteStart and WriteEnd
            dataflow.pipepool.(nextnode)(1).WriteStart = dataflow.pipepool.(nextnode)(1).WritePoint;
            dataflow.pipepool.(nextnode)(1).WriteEnd = dataflow.pipepool.(nextnode)(1).WriteStart + Nviewcurrshot - 1;
            % reset nextpool's ReadPoint
            dataflow.pipepool.(nextnode)(1).ReadPoint = dataflow.pipepool.(nextnode)(1).WritePoint;
            dataflow.pipepool.(nextnode)(1).ReadStart = dataflow.pipepool.(nextnode)(1).ReadPoint;
            dataflow.pipepool.(nextnode)(1).ReadEnd = dataflow.pipepool.(nextnode)(1).ReadStart + Nviewcurrshot - 1;
            % reset nextpool's AvialPoint
            dataflow.pipepool.(nextnode)(1).AvailPoint = dataflow.pipepool.(nextnode)(1).ReadPoint-1;
%             % close the shotstart
%             dataflow.pipepool.(nextnode)(1).isshotstart = false;
            % check poolsize
            if dataflow.pipepool.(nextnode)(1).circulatemode
                if dataflow.pipepool.(nextnode)(1).poolsize ~= Nviewcurrshot
                    % should be some buffer re-malloc here.
                    dataflow.pipepool.(nextnode)(1).poolsize = Nviewcurrshot; 
                end
            end
            % Here is a check point for C++ codes to new or re-malloc the buffers.
            
            % ishot++ 
            dataflow.buffer.(nodename).ishot = dataflow.buffer.(nodename).ishot + 1;
        end
        % nextpoolleft
        nextpoolleft = poolpspaceleft(dataflow.pipepool.(nextnode)(1));
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
    % set torunpoststep = true
    status.currentjob.torunpoststep = true;

    if writenum == 0 && AvailNumber>0
        % to pass
        status.jobdone = 6;
    elseif writenum < AvailNumber
        status.jobdone = 2;
    end
    % I know, else while writenum=AvailNumber>0, status.jobdone=1 if iblock > Nblock, or status.jobdone=2

    % return while to pass or error
    if status.jobdone == 6 || status.jobdone == 0
        return;
    end

    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);

end
% done

end


function [currpool, currdata, offset, headprm, outcfg] = loadrawdata(currpool, currdata, filename, ...
    prmflow, startview, viewnum, circulatereading)
% load rawdata (or offset) from the file filename to dataflow (or currdata)
% Note that not all the tags in the rawdata can be read to the dataflow, see function raw2dataflow, only a limited part of them
% will be extracted to support most cali and recon tasks.

if nargin < 7
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
    [~, fileNAME, fileEXT] = fileparts(filename);
end

switch lower(fileEXT)
    case {'.raw', '.bin'}
        % It is like, raw = loaddata(filename, IOstandard);
        % But not enough, yep, let'd do it
        rawcfg = prmflow.raw.formatcfg;
        versionflag = prmflow.raw.formatversion;
        if isempty(rawcfg)
            [cfgfile, versionflag] = cfgmatchrule(filename, IOstandard);
            rawcfg = readcfgfile(cfgfile);
        end
        rawcfg.number = viewnum;
        viewskip = startview-1;
        if circulatereading && isfield(rawcfg, 'filelength')
            viewskip = mod(viewskip, rawcfg.filelength);
        end
        fid = fopen(filename, 'r');
        [raw, outcfg] = sparsepack(fid, rawcfg, viewskip);
        readingnumber = outcfg.number;
        % circulate reading
        if circulatereading && ~isinf(viewnum)
            while readingnumber < viewnum
                outcfg.number = viewnum - readingnumber;
                fseek(fid, 0, 'bof');
                [raw1, outcfg] = sparsepack(fid, outcfg, 0);
                raw = cat(2, raw, raw1);
                readingnumber = readingnumber + outcfg.number;
            end
        end
        fclose(fid);
        if isinf(viewnum)
            viewnum = readingnumber;
        elseif viewnum > readingnumber
            % unexpected file end
            1;
            % should have a warning
        end

        % version flag
        keyversion = regexp(versionflag(2:end), '\.', 'split');
        if length(keyversion) > 1
            keyversion = str2double(keyversion{end-1});
        elseif ~isempty(keyversion{1})
            keyversion = str2double(keyversion{1});
        else
            keyversion = 0;
        end

        % to dataflow
        if keyversion == 2
            % external rawdata sparse for deployments
            [rawhead, rawdata, headprm] = raw2dataflowV2(raw);
        else
            % v1.* or v0.* (for general using)
            [rawhead, rawdata, headprm] = raw2dataflow(raw);
        end
        % save the versionflag
        headprm.versionflag = versionflag;
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
rawhead = defaultrawhead(rawhead, prmflow, startview, viewnum);

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
rawhead.AngleEncoder = [raw.AngleEncoder];
rawhead.ReadingNumber = [raw.ReadingNumber];
rawhead.IntegrationTime = [raw.IntegrationTime];
rawhead.ShotNumber = [raw.ShotNumber];
% rawhead.TimeStamp = [raw.TimeStamp];
rawhead.mA = single([raw.mA]);
rawhead.KV = single([raw.KV]);
rawhead.TableEncoder = [raw.TableEncoder];
if isfield(raw, 'TableGear')
    rawhead.TableGear = [raw.TableGear];
end
% We may have a list to configure in reading those tags to rawhead, which could be depending on protocol.

% rawdata
rawdata = single([raw.RawData]);
% I know all the fields of rawhead and rawdata are in size n*Nview

if ~isempty(raw)
    headprm.PackageVersion = raw(1).PackageVersion;
    headprm.SeriesNumber = raw(1).SeriesNumber;
    headprm.StartSlice = raw(1).StartSlice;
    headprm.EndSlice = raw(1).EndSlice;
    headprm.SliceMergescale = raw(1).SliceMergescale;
    headprm.SliceNumber = raw(1).SliceNumber;
    headprm.RawDataSize = raw(1).RawDataSize;
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

function rawhead = defaultrawhead(rawhead, prmflow, startview, viewnum)
% hard coded default rawhead
% to fill up the missing fields and align the data class

rawheadfields = prmflow.raw.rawheadfields;
angulationcode = prmflow.system.angulationcode;
angulationzero = prmflow.system.angulationzero;
viewpershot = prmflow.raw.viewpershot;
couchspeed = prmflow.raw.couchspeed;

% viewindex and shotindex
viewindex = double(startview) : double(startview+viewnum)-1;
viewpershotCUM = [0; cumsum(viewpershot(:))];
shotindex = sum(viewindex > viewpershotCUM, 1);

% to fill up the ReadingNumber if not exist
if ~isfield(rawhead, 'ReadingNumber') && isincell(rawheadfields, 'ReadingNumber')
    rawhead.ReadingNumber = uint32(viewindex - viewpershotCUM(shotindex)');
end

% to fill up the ShotNumber if not exist
if ~isfield(rawhead, 'ShotNumber') && isincell(rawheadfields, 'ShotNumber')
    rawhead.ShotNumber = uint16(shotindex);
end
% to fill up the ShotStart if not exist
if ~isfield(rawhead, 'ShotStart') && isincell(rawheadfields, 'ShotStart')
    currShotStart = viewindex - viewpershotCUM;
    rawhead.ShotStart = int16(any(currShotStart==1, 1)) - int16(any(currShotStart==0, 1));
end
% the ShotStart is demanded to be 1 for the first view of a shot, -1 for the end view of a shot and 0 for the others. 
% Note: in pipeline console we will not use the rawhead.ShotStart to judge the shot start/end. Actually we never use any
% information from the the dataflow to govern the pipeline accesses.

% refblock (should be asked by referencecorr)
if ~isfield(rawhead, 'refblock') && isincell(rawheadfields, 'refblock')
    rawhead.refblock = false(2, viewnum);
end

% TableGear (should be asked by Backprojection for Conveyor)
if ~isfield(rawhead, 'TableGear') && isincell(rawheadfields, 'TableGear')
    rawhead.TableGear = int16(ones(1, viewnum).*sign(couchspeed));
end

% viewangle (should be asked by some calibration nodes or never)
if ~isfield(rawhead, 'viewangle') && isincell(rawheadfields, 'viewangle')
    rawhead.viewangle = single((double(rawhead.AngleEncoder) - angulationzero)./angulationcode.*(pi*2));
end

% remove/fill the fields
for field_ii = fieldnames(rawhead)'
    if ~isincell(rawheadfields, field_ii{1})
        rawhead = rmfield(rawhead, field_ii{1});
    end
end
for field_ii = rawheadfields
    if ~isfield(rawhead, field_ii)
        rawhead.(field_ii{1}) = zeros(1, viewnum);
    end
end

% data class
if isincell(rawheadfields, 'ReadingNumber') && ~isa(rawhead.ReadingNumber, 'uint32')
    rawhead.ReadingNumber = uint32(rawhead.ReadingNumber);
end
if isincell(rawheadfields, 'ShotNumber') && ~isa(rawhead.ShotNumber, 'uint16')
    rawhead.ShotNumber = uint16(rawhead.ShotNumber);
end
if isincell(rawheadfields, 'ShotStart') && ~isa(rawhead.ShotStart, 'int16')
    rawhead.ShotStart = int16(rawhead.ShotStart);
end
if isincell(rawheadfields, 'AngleEncoder') && ~(isa(rawhead.AngleEncoder, 'uint32') || isa(rawhead.AngleEncoder, 'uint64'))
    rawhead.AngleEncoder = uint32(rawhead.AngleEncoder);
end
if isincell(rawheadfields, 'TableEncoder') && ~(isa(rawhead.TableEncoder, 'uint32') || isa(rawhead.TableEncoder, 'uint64'))
    rawhead.TableEncoder = uint32(rawhead.TableEncoder);
end
if isincell(rawheadfields, 'TableGear') && ~isa(rawhead.TableGear, 'int16')
    rawhead.TableGear = int16(rawhead.TableGear);
end
if isincell(rawheadfields, 'viewangle') && ~isa(rawhead.viewangle, 'single')
    rawhead.viewangle = single(rawhead.viewangle);
end
if isincell(rawheadfields, 'refblock') && ~isa(rawhead.refblock, 'logical')
    rawhead.refblock = logical(rawhead.refblock);
end
if isincell(rawheadfields, 'moverdis') && ~isa(rawhead.viewangle, 'double')
    rawhead.moverdis = double(rawhead.moverdis);
end

% other
% to fill up viewangle by AngleEncoder (We will remove this soon), do not use it on deployment
if isfield(rawhead, 'AngleEncoder')
    % the rawhead.viewangle is still used in calibration nodes, though we omited it in recon nodes.
    rawhead.viewangle = single((double(rawhead.AngleEncoder) - angulationzero)./angulationcode.*(pi*2));
    % we will put the asks in the calibration nodes
end

end

function [buffer, dataflow, prmflow, status] = fakekernelreadrawdata(buffer, dataflow, prmflow, status)
% fake read rawdata

% nodename
nodename = status.nodename;
% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% shots to read
shotnum = prmflow.raw.Nshot;
startshot = prmflow.raw.startshot;
viewpershot = prmflow.raw.viewpershot;
if isfield(buffer, 'datablock_onoff') && buffer.datablock_onoff
    % to read data by blocks
    datablocksize = buffer.datablocksize;
    iblock = buffer.iblock;
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

% fake rawdata
if ~isfield(dataflow, 'rawdata')
    dataflow.rawdata = [];
end
% fake rawhead
if ~isfield(dataflow, 'rawhead')
    dataflow.rawhead = struct();
end
if ~isfield(dataflow.rawhead, 'Reading_Number')
    dataflow.rawhead.Reading_Number = [];
end
if ~isfield(dataflow.rawhead, 'Shot_Number')
    dataflow.rawhead.Shot_Number = [];
end
currReading_Number = startview : startview+viewnum-1;
currShotnumber = ceil(currReading_Number / viewpershot);

if pipeline_onoff
    writeindex = buffer.outputpool.WritePoint : buffer.outputpool.WritePoint + viewnum - 1;
    dataflow.rawdata(:, writeindex) = zeros(1, viewnum, 'single');
    dataflow.rawhead.Reading_Number(:, writeindex) = currReading_Number;
    dataflow.rawhead.Shot_Number(:, writeindex) = currShotnumber;
else
    dataflow.rawdata = [dataflow.rawdata, zeros(1, viewnum, 'single')];
    dataflow.rawhead.Reading_Number = [dataflow.rawhead.Reading_Number currReading_Number];
    dataflow.rawhead.Shot_Number = [dataflow.rawhead.Shot_Number currShotnumber];
end

% datablock
if isfield(buffer, 'datablock_onoff') && buffer.datablock_onoff
    % contol the looping of data blocks
    iblock = buffer.iblock;
    Nblock = buffer.Nblock;
    if iblock < Nblock
        status.jobdone = 2;
    else
%         status.pipeline.(nodename).sleeping = true;
        status.jobdone = 1;
    end
    % iblock++
    buffer.iblock = iblock + 1;
    % WritePoint for the 'hidden' node
    buffer.outputpool.WritePoint = buffer.outputpool.WritePoint + viewnum;
    buffer.outputpool.AvailNumber = buffer.outputpool.AvailNumber + viewnum;
else
    status.jobdone = true;
end

if pipeline_onoff && ~buffer.datablock_onoff
    % What? Shouldn't we forbid this behavior?
    buffer.outputpool.WritePoint = buffer.outputpool.WritePoint + viewnum;
    buffer.outputpool.AvailNumber = buffer.outputpool.AvailNumber + viewnum;
end

end



function [dataflow, prmflow, status] = fakenodereadrawdata(dataflow, prmflow, status)

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% % kill the output pool
% if isfield(dataflow.buffer.(nodename), 'outputpool')
%     dataflow.buffer.(nodename) = rmfield(dataflow.buffer.(nodename), 'outputpool');
% end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% shots to read
shotnum = prmflow.raw.Nshot;
startshot = prmflow.raw.startshot;
viewpershot = prmflow.raw.viewpershot;
% if prmflow.raw.datablock_onoff
% to read data by blocks
datablocksize = prmflow.raw.datablocksize;
iblock = prmflow.raw.iblock;
Nblock = prmflow.raw.Nblock;
if iblock <= Nblock
    readingdata = true;
    if length(viewpershot) == 1
        startview = sum(datablocksize(1:iblock-1)) + (startshot-1)*viewpershot + 1;
    else
        startview = sum(datablocksize(1:iblock-1)) + sum(viewpershot(1 : startshot-1)) + 1;
    end
    viewnum = datablocksize(iblock);

    % iblock++
    prmflow.raw.iblock = iblock + 1;
else
    % to pass
    readingdata = false;
end
% else
%     % read all data in one block
%     if length(viewpershot) == 1
%         startview = (startshot-1)*viewpershot + 1;
%         viewnum = shotnum * viewpershot;
%     else
%         startview = sum(viewpershot(1 : startshot-1)) + 1;
%         viewnum = sum(viewpershot(startshot : startshot+shotnum-1));
%     end
% end

% fake rawdata
if readingdata
    if pipeline_onoff
        [dataflow.pipepool.(nodename), dataflow.pipepool.(nodename).data] = ...
            readfakeraw(dataflow.pipepool.(nodename), dataflow.pipepool.(nodename).data, startview, viewnum, viewpershot);
    else
        [~, dataflow] = readfakeraw([], dataflow, startview, viewnum, viewpershot);
    end
end

% datablock
% if prmflow.raw.datablock_onoff
    % contol the looping of data blocks   
    if prmflow.raw.iblock <= Nblock
        % keep waking
        status.jobdone = 2;
    else
%         status.pipeline.(nodename).sleeping = true;
        status.jobdone = 1;
    end
    
    % WritePoint for the 'hidden' node
% else
%     status.jobdone = true;
% end
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
            % view number (Nend) of current shot
            ishot = dataflow.buffer.(nodename).ishot + startshot;
            if length(viewpershot) == 1
                Nend = viewpershot;
            else
                Nend = viewpershot(ishot);
            end
            % reset currpool's ReadStart and ReadEnd
            dataflow.pipepool.(nodename).ReadStart = dataflow.pipepool.(nodename).ReadPoint;
            dataflow.pipepool.(nodename).ReadEnd = dataflow.pipepool.(nodename).ReadStart + Nend - 1;
            % reset nextpool's WriteStart and WriteEnd
            dataflow.pipepool.(nextnode).WriteStart = dataflow.pipepool.(nextnode).WritePoint;
            dataflow.pipepool.(nextnode).WriteEnd = dataflow.pipepool.(nextnode).WriteStart + Nend - 1;
            % reset nextpool's ReadPoint
            dataflow.pipepool.(nextnode).ReadPoint = dataflow.pipepool.(nextnode).WritePoint;
            dataflow.pipepool.(nextnode).ReadStart = dataflow.pipepool.(nextnode).ReadPoint;
            dataflow.pipepool.(nextnode).ReadEnd = dataflow.pipepool.(nextnode).ReadStart + Nend - 1;
            % reset nextpool's AvialPoint
            dataflow.pipepool.(nextnode).AvailPoint = dataflow.pipepool.(nextnode).ReadPoint-1;
            % check poolsize
            if dataflow.pipepool.(nextnode).circulatemode
                if dataflow.pipepool.(nextnode).poolsize ~= viewpershot
                    % should be some buffer re-malloc here.
                    dataflow.pipepool.(nextnode).poolsize = viewpershot; 
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

        % copy rawdata to next pool
        writenum = min(nextpoolleft, AvailNumber);
        dataflow.pipepool.(nextnode).data = pooldatacopy(dataflow.pipepool.(nodename), dataflow.pipepool.(nodename).data, ...
            dataflow.pipepool.(nextnode), dataflow.pipepool.(nextnode).data, writenum, [], true);
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
    
    if writenum == 0 && AvailNumber>0
        status.jobdone = 6;
    elseif writenum < AvailNumber
        status.jobdone = 2;
    end
    % I know, else while writenum=AvailNumber>0, status.jobdone=1 if iblock > Nblock, or status.jobdone=2

%     % dataflow mirror (useless)
%     if isfield(nodeprm, 'dataflowmirror') && nodeprm.dataflowmirror
%         dataflow.rawdata = dataflow.pipepool.(nodename).rawdata;
%         dataflow.rawhead = dataflow.pipepool.(nodename).rawhead;
%     end
end
% done

status.errorcode = 0;
status.errormsg = [];

end


function [currpool, currdata] = readfakeraw(currpool, currdata, startview, viewnum, viewpershot)
% to 'read' the fake data.
% in non-pipeline mode, the currdata shall be dataflow and currpool shall be empty or a fake pipepool in
% dataflow.buffer.(nodename)

% fake rawdata
if ~isfield(currdata, 'rawdata')
    currdata.rawdata = [];
end
% fake rawhead
if ~isfield(currdata, 'rawhead')
    currdata.rawhead = struct();
end
if ~isfield(currdata.rawhead, 'Reading_Number')
    currdata.rawhead.Reading_Number = [];
end
if ~isfield(currdata.rawhead, 'Shot_Number')
    currdata.rawhead.Shot_Number = [];
end
% Reading_Number and Shotnumber
viewindex = startview : startview+viewnum-1;
currReading_Number = mod(viewindex-1, viewpershot(1))+1;
currShotnumber = ceil(viewindex / viewpershot(1));  % not so strict
% ShotStart
currShotStart = mod(viewindex, viewpershot(1));
currShotStart(currShotStart==0) = -1;
currShotStart(currShotStart>1) = 0;
currShotStart = int16(currShotStart);

if ~isempty(currpool)
    writeindex = currpool.WritePoint : currpool.WritePoint + viewnum - 1;
    currdata.rawdata(:, writeindex) = zeros(1, viewnum, 'single') + 1;
    currdata.rawhead.Reading_Number(:, writeindex) = currReading_Number;
    currdata.rawhead.Shot_Number(:, writeindex) = currShotnumber;
    currdata.rawhead.Shot_Start(:, writeindex) = currShotStart;
    
    currpool.WritePoint = currpool.WritePoint + viewnum;
    currpool.AvailPoint = currpool.AvailPoint + viewnum;
else
    currdata.rawdata = [currdata.rawdata, zeros(1, viewnum, 'single')];
    currdata.rawhead.Reading_Number = [currdata.rawhead.Reading_Number currReading_Number];
    currdata.rawhead.Shot_Number = [currdata.rawhead.Shot_Number currShotnumber];
    currdata.rawhead.Shot_Start = [currdata.rawhead.ShotStart currShotStart];
end

currdata.rawhead.OrigReadingNumber = currdata.rawhead.Reading_Number;
end

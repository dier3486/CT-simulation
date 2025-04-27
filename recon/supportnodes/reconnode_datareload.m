function [dataflow, prmflow, status] = reconnode_datareload(dataflow, prmflow, status)
% support node, to reload data from dataflow to next pool while pipeline-on
%   [dataflow, prmflow, status] = reconnode_datareload(dataflow, prmflow, status);
% So called 'reload' is to copy data by blocks from dataflow (but not files) to next node's input pool to reload the pipeline.
% 

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);


if ~status.pipeline.(nodename).pipeline_onoff
    warning('reconnode_datareload can only work in pipeline mode!');
    return;
    % in non-pipeline mode the nodes can visit the data in dataflow directly that needs not 'reload'.
end

% active (in working)
if isfield(status.pipeline.(nodename), 'active') && ~status.pipeline.(nodename).active
    % pass
    status.jobdone = 3;
    return;
end

reloadfields = nodeprm.reloadfields;
reloadfields = reloadfields(isfield(dataflow, reloadfields));
% We shall not copy all fields in dataflow (no 'all' option), because really all data are in dataflow, all the private buffer
% all the pools, all the backups and all the outputs.

% ini pipepool
if dataflow.buffer.(nodename).iblock == 1
    % rawlength
    if any(strcmp(reloadfields, 'rawdata'))
        rawlength = size(dataflow.rawdata, 2);
    elseif any(strcmp(reloadfields, 'rawhead'))
        headfields = fieldnames(dataflow.rawhead);
        rawlength = size(dataflow.rawhead.(headfields{1}), 2);
    else
        rawlength = size(dataflow.(reloadfields{1}), 2);
    end
    % check error
    if isinf(dataflow.buffer.(nodename).Nview)
        dataflow.buffer.(nodename).Nview = rawlength;
    elseif dataflow.buffer.(nodename).Nview ~= rawlength
        % error
        error('View number unmatch!');
    end
    if isinf(dataflow.buffer.(nodename).viewpershot)
        dataflow.buffer.(nodename).viewpershot = rawlength;
        if dataflow.buffer.(nodename).Nshot ~= 1
            error('Missed shotend!');
        end
    end
    % ini ReadEnd
    dataflow.pipepool.(nodename).WriteEnd = dataflow.pipepool.(nodename).ReadPoint + rawlength - 1;
    dataflow.buffer.(nodename).Nblock = ceil(rawlength / nodeprm.datablock);
    % reset poolsize
    dataflow.pipepool.(nodename).poolsize = rawlength;
end
1;
% move a block
dataflow.pipepool.(nodename).AvailPoint = min(dataflow.pipepool.(nodename).AvailPoint + nodeprm.datablock, ...
    dataflow.pipepool.(nodename).WriteEnd);
% iblock++
if dataflow.buffer.(nodename).iblock < dataflow.buffer.(nodename).Nblock
    dataflow.buffer.(nodename).iblock = dataflow.buffer.(nodename).iblock + 1;
    status.jobdone = 2;
else
    status.jobdone = 1;
end

% view number per shot
viewpershot = dataflow.buffer.(nodename).viewpershot;

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
        ishot = dataflow.buffer.(nodename).ishot + 1;
        if length(viewpershot) == 1
            Nviewcurrshot = viewpershot;
        else
            Nviewcurrshot = viewpershot(ishot);
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

    % copy rawdata to next pool
    writenum = min(nextpoolleft, AvailNumber);
    dataflow.pipepool.(nextnode).data = pooldatacopy(dataflow.pipepool.(nodename), dataflow, ...
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
% else, see line 64 % iblock++

% post step
[dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);

% clear data
if status.jobdone==1 && nodeprm.clearafterreload
    dataflow = rmfield(dataflow, reloadfields);
end

status.errorcode = 0;
status.errormsg = [];
end
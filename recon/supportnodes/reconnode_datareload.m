function [dataflow, prmflow, status] = reconnode_datareload(dataflow, prmflow, status)
% support node, to reload data from dataflow to next pool
%   [dataflow, prmflow, status] = reconnode_datareload(dataflow, prmflow, status);
% So called 'reload' is to copy data in blocks from dataflow (but not files) to next node's input pool.
% 

if ~status.pipeline_onoff
    warning('reconnode_datareload can only work in pipeline mode!');
    return;
    % in non-pipeline mode the nodes can visit the data in dataflow directly that needs not any 'reload'.
end

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);
nextnode = status.pipeline.(nodename).nextnode;

% active (in working)
if isfield(status.pipeline.(nodename), 'active')
    active = status.pipeline.(nodename).active;
else
    active = true;
end

% datablock (size)
if isfield(nodeprm, 'datablock')
    datablock = nodeprm.datablock;
else
    datablock = 512;
end
% data fields
if isfield(nodeprm, 'reloadfields')
    reloadfields = regexp(nodeprm.reloadfields, '(, +)|(,)', 'split');
else
    % default reloadfields are rawdata and rawhead
    reloadfields = {'rawdata', 'rawhead'};
end
reloadfields = reloadfields(isfield(dataflow, reloadfields));
% We shall not copy all fields in dataflow (no 'all' option), because really all data are in dataflow, all the private buffer
% all the pools, all the backups and all the outputs.

% datalength (mostly =Nview)
if any(strcmp(reloadfields, 'rawdata'))
    rawlength = size(dataflow.rawdata, 2);
elseif any(strcmp(reloadfields, 'rawhead'))
    headfields = fieldnames(dataflow.rawhead);
    rawlength = size(dataflow.rawhead.(headfields{1}), 2);
else
    rawlength = size(dataflow.(reloadfields{1}), 2);
end

% ini iblock
if ~isfield(dataflow.buffer.(nodename), 'iblock')
    dataflow.buffer.(nodename).iblock = 1;
end
% ini pipepool
if ~isfield(dataflow.buffer.(nodename), 'pipepool')
    dataflow.buffer.(nodename).pipepool.ReadPoint = 1;
    dataflow.buffer.(nodename).pipepool.WritePoint = rawlength + 1;
end

% the current pool of the datareload node is not in status.pipepool, but in private buffer
currpool = dataflow.buffer.(nodename).pipepool;

% copy data
if ~isempty(nextnode)
    if active
        statusnext = status.pipepool.(nextnode);
        % data size to be written in the next pool
        datasize = min(datablock, currpool.WritePoint - currpool.ReadPoint);
        writenum = min(datasize, min(statusnext.WriteEnd, statusnext.poolsize) - statusnext.WritePoint + 1);
        1;
        % copy data
        [dataflow.pipepool.(nextnode), writenum] = pooldatacopy(dataflow, dataflow.pipepool.(nextnode), ...
            currpool.ReadPoint, statusnext.WritePoint, writenum, reloadfields, true);

        % move next pool's write point
        status.pipepool.(nextnode).WritePoint = status.pipepool.(nextnode).WritePoint + writenum;
        % move current pool's read point
        dataflow.buffer.(nodename).pipepool.ReadPoint = dataflow.buffer.(nodename).pipepool.ReadPoint + writenum;

        % recycle
        % not yet

        % job done
        if dataflow.buffer.(nodename).pipepool.ReadPoint < dataflow.buffer.(nodename).pipepool.WritePoint
            % technically done
            status.jobdone = 2;
        else
            % job done, sleep
            status.jobdone = 1;
        end
    else % in working..
        % do nothing
        [dataflow, prmflow, status] = reconnode_donothing(dataflow, prmflow, status);
    end
else
    % move current pool's read point
    datasize = min(datablock, currpool.WritePoint - currpool.ReadPoint);
    dataflow.buffer.(nodename).pipepool.ReadPoint = dataflow.buffer.(nodename).pipepool.ReadPoint + datasize;
    % job done
    if dataflow.buffer.(nodename).pipepool.ReadPoint < dataflow.buffer.(nodename).pipepool.WritePoint
        % technically done
        status.jobdone = 2;
    else
        % job done, sleep
        status.jobdone = 1;
    end
end

status.errorcode = 0;
status.errormsg = [];
end
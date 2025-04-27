function [dataflow, prmflow, status] = reconnode_databuffer(dataflow, prmflow, status)
% support node, to buffer data
% [dataflow, prmflow, status] = reconnode_databuffer(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;
if ~pipeline_onoff
    error('reconnode_databuffer can only run in pipeline mode!');
end

% the data to be buffered
if isfield(nodeprm, 'alldata')  && nodeprm.alldata
    bufferfields = {};
    % {} will buffer all fields, (and {''} will buffer nothing)
    if isfield(nodeprm, 'excludefields')
        % TBC
        1;
    end
else
    if isfield(nodeprm, 'bufferfields')
        bufferfields = regexp(nodeprm.bufferfields, '(, +)|(,)', 'split');
    else
        bufferfields = {''};
        % if that, nothing will be buffered
    end
end
% TBC

% if to copy the buffer back to dataflow (after previous works done)
if isfield(nodeprm, 'copytodataflow')
    copytodataflow = nodeprm.copytodataflow;
else
    copytodataflow = false;
end
% if to stuck the pipeline until previous works done, used to link a pipeline node with a non-pipeline node.
if isfield(nodeprm, 'stuck')
    stuck_flag = nodeprm.stuck;
else
    stuck_flag = false;
end

% these should be initialed (we may have a reconnode_databufferprepare.m)
if ~isfield(dataflow.buffer.(nodename), 'data')
    dataflow.buffer.(nodename) = status.defaultpool;
    dataflow.buffer.(nodename).datafields = dataflow.pipepool.(nodename).datafields;
    dataflow.buffer.(nodename).data = struct();
    for ifield = dataflow.pipepool.(nodename).datafields
        if isstruct(dataflow.pipepool.(nodename).data.(ifield{1}))
            dataflow.buffer.(nodename).data.(ifield{1}) = struct();
        else
            dataflow.buffer.(nodename).data.(ifield{1}) = dataflow.pipepool.(nodename).data.(ifield{1})(:, []);
        end
    end
    1;
end
% if ~isfield(dataflow.buffer.(nodename), 'ReadPoint')
%     dataflow.buffer.(nodename).ReadPoint = 1;
% end
% if ~isfield(dataflow.buffer.(nodename), 'WritePoint')
%     dataflow.buffer.(nodename).WritePoint = 1;
% end

% copy data from currpool to buffer
% [status.pipepool.(nodename), dataflow.buffer.(nodename), dataflow.buffer.(nodename).data, Rinfo] = ...
%     pooldataoutput_cyc2(status.pipepool.(nodename), dataflow.pipepool.(nodename), ...
%     dataflow.buffer.(nodename), dataflow.buffer.(nodename).data);
[dataflow.pipepool.(nodename), dataflow.buffer.(nodename), Rinfo] = ...
    pooldataoutput(dataflow.pipepool.(nodename), dataflow.buffer.(nodename));
% add Avail
dataflow.buffer.(nodename).AvailPoint = dataflow.buffer.(nodename).AvailPoint + Rinfo.writenum;

% WritePoint = dataflow.buffer.(nodename).WritePoint;
% writenum = max(0, statuspool.WritePoint - statuspool.ReadPoint);
% % copy data from pipepool.(nodename) to buffer.(nodename).data
% [dataflow.buffer.(nodename).data, writenum] = pooldatacopy(dataflow.pipepool.(nodename), ...
%     dataflow.buffer.(nodename).data, statuspool.ReadPoint, WritePoint, writenum, bufferfields, true);
% % move buffer WritePoint
% dataflow.buffer.(nodename).WritePoint = dataflow.buffer.(nodename).WritePoint + writenum;
% % move pool ReadPoint
% if stuck_flag
%     status.pipepool.(nodename).ReadPoint = status.pipepool.(nodename).ReadPoint + writenum;
% end
% % else the ReadPoint will be moved by reconnode_donothing

if ~stuck_flag
    % copy the buffer to nextpool
    [dataflow.buffer.(nodename), dataflow.pipepool.(nextnode), Rinfo] = ...
        pooldataoutput(dataflow.buffer.(nodename), dataflow.pipepool.(nextnode));
    % jobdone
    status.jobdone = Rinfo.jobflag;
end

% all the previous node sleeping?
nodesindex = 0 : status.pipeline.(nodename).index - 1;
allsleeping = nodestatecheck(status.pipeline, 'sleeping', nodesindex, 'all');
% if to copy buffer data to dataflow
if copytodataflow && allsleeping
    % copy all the buffer data to dataflow (overwrite!)
    dataflow.buffer.(nodename).ReadPoint = 1;
    writenum = dataflow.buffer.(nodename).WritePoint - 1;
    [dataflow, ~] = pooldatacopy(dataflow.buffer.(nodename), dataflow.buffer.(nodename).data, ...
        status.defaultpool, dataflow, writenum, [], true);

    % clear buffer
    dataflow.buffer.(nodename).data = poolclear(dataflow.buffer.(nodename).data);
    dataflow.buffer.(nodename).ReadPoint = 1;
    dataflow.buffer.(nodename).WritePoint = 1;
    dataflow.buffer.(nodename).AvailPoint = 0;
%     % series done
%     status.pipeline.(nodename).seriesdone = true;
end

if stuck_flag
    % job done
    if allsleeping
        status.jobdone = 1;
        % 1: to wake up the next node which is supposed a non-pipeline node.
    else
        if Rinfo.writenum > 0
            status.jobdone = 4;
            % 4: technically pass, don't wake up the next node.
        else
            status.jobdone = 3;
            % 3: pass, did nothing
        end
    end
end

status.errorcode = 0;
status.errormsg = [];
end
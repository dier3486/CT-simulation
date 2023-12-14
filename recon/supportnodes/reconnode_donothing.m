function [dataflow, prmflow, status] = reconnode_donothing(dataflow, prmflow, status)
% support node, do nothing
% [dataflow, prmflow, status] = reconnode_donothing(dataflow, prmflow, status);
% It is the most important node on pipeline ;) which is a template of the pipeline nodes.

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

if pipeline_onoff
    % do nothing is not trivial in pipeline
    nextnode = status.pipeline.(nodename).nextnode;
    if ~isempty(nextnode)
        % to copy the data in the pipepool to next node
        statuscurr = status.pipepool.(nodename);
        statusnext = status.pipepool.(nextnode);
        datasize = max(0, statuscurr.WritePoint - statuscurr.ReadPoint);
        writenum = min(datasize, min(statusnext.WriteEnd, statusnext.poolsize) - statusnext.WritePoint + 1);
        % call pooldatacopy
        [dataflow.pipepool.(nextnode), writenum] = pooldatacopy(dataflow.pipepool.(nodename), ...
            dataflow.pipepool.(nextnode), statuscurr.ReadPoint, statusnext.WritePoint, writenum);
        % move next pool's write point
        status.pipepool.(nextnode).WritePoint = status.pipepool.(nextnode).WritePoint + writenum;
        % move current pool's read point
        status.pipepool.(nodename).ReadPoint = status.pipepool.(nodename).ReadPoint + writenum;
    else
        % Even no next node to output, we also need to move the ReadPoint in the current pool to help consol node releasing the
        % buffer resources.
        statuscurr = status.pipepool.(nodename);
        datasize = max(0, statuscurr.WritePoint - statuscurr.ReadPoint);
        writenum = datasize;
        % move current pool's read point
        status.pipepool.(nodename).ReadPoint = status.pipepool.(nodename).ReadPoint + writenum;
    end
    % that is a sandard data transfer mode for pipeline nodes.
    if datasize == 0
        % donothing did nothing, pass
        status.jobdone = 3;
    elseif writenum < datasize
        % done and keep waking
        status.jobdone = 2;
    else
        % normally done
        status.jobdone = 1;
    end
else
    % really do nothing
    status.jobdone = true;
end

status.errorcode = 0;
status.errormsg = [];
end
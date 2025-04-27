function [dataflow, prmflow, status] = reconnode_donothing(dataflow, prmflow, status)
% support node, do nothing
% [dataflow, prmflow, status] = reconnode_donothing(dataflow, prmflow, status);
% It is the most important node on pipeline ;) which is a template of the pipeline nodes.

% parameters set in pipe
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% prio step
if pipeline_onoff
    % do nothing is not trivial in pipeline (it is)
    
    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    
    if status.currentjob.topass
        % error or pass
        return;
    end
else
    status.jobdone = true;
end

% node kernel function
1;
% nothing

% post step
if pipeline_onoff
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end

status.errorcode = 0;
status.errormsg = [];
end

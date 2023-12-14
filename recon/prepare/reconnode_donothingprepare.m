function [dataflow, prmflow, status] = reconnode_donothingprepare(dataflow, prmflow, status)
% prepare node, prepare to do nothing
% [dataflow, prmflow, status] = reconnode_donothingprepare(dataflow, prmflow, status);
% wow~

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% input pool
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpool;
    % We shall allow user to define the input pool fields
end

status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end
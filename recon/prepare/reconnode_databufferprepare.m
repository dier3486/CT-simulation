function [dataflow, prmflow, status] = reconnode_databufferprepare(dataflow, prmflow, status)
% prepare node, data buffer prepare
% [dataflow, prmflow, status] = reconnode_databufferprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% pipeline_onoff
if isfield(nodeprm, 'pipeline_onoff')
    pipeline_onoff = status.pipeline_onoff & nodeprm.pipeline_onoff;
else
    pipeline_onoff = status.pipeline_onoff;
end

% input pool
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpooldata;
    % We shall allow user to define the input pool fields
    dataflow.buffer.(nodename) = struct();
end

status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end
function [dataflow, prmflow, status] = reconnode_pipelinestuckprepare(dataflow, prmflow, status)
% prepare node, pipeline stuck prepare
% [dataflow, prmflow, status] = reconnode_pipelinestuckprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% pipeline stuck is a hard coded reconnode_databuffer
if ~isfield(nodeprm, 'alldata') && ~isfield(nodeprm, 'bufferfields')
    prmflow.pipe.(nodename).alldata = true;
end
% force
prmflow.pipe.(nodename).copytodataflow = true;
prmflow.pipe.(nodename).stuck = true;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% input pool
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpool;
    % We shall allow user to define the input pool fields
    dataflow.buffer.(nodename) = struct();
end

status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end
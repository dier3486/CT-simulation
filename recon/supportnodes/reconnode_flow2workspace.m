function [dataflow, prmflow, status] = reconnode_flow2workspace(dataflow, prmflow, status)
% virus node, output dataflow, prmflow and status to workspace
% [dataflow, prmflow, status] = reconnode_flow2workspace(dataflow, prmflow, status);


% parameters set in pipe
nodeprm = prmflow.pipe.(status.nodename);

% backup index
if isfield(nodeprm, 'index')
    flowindex = num2str(nodeprm.index);
else
    flowindex = num2str(status.seriesindex);
end

if ~isempty(flowindex)
    flowindex = ['_' flowindex];
end

% dataflow, prmflow and status
assignin('base', ['dataflow' flowindex], dataflow);
assignin('base', ['prmflow' flowindex], prmflow);
assignin('base', ['status' flowindex], status);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end

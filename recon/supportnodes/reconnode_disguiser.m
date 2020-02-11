function [dataflow, prmflow, status] = reconnode_disguiser(dataflow, prmflow, status)
% virus node, I am your mirror
% [dataflow, prmflow, status] = reconnode_disguiser(dataflow, prmflow, status);
% WARN: uncontrollable

% parameters set in pipe
disgprm = prmflow.pipe.(status.nodename);

% I know the status.nodename is 'disguiser' or 'disguiser_*'
if isfield(disgprm, 'node') && ~isempty(disgprm.node)
    myfun = str2func(['reconnode_' disgprm.node]);
    [dataflow, prmflow, status] = myfun(dataflow, prmflow, status);
else
    % do nothing
    % status
    status.jobdone = true;
    status.errorcode = 0;
    status.errormsg = [];
end

end
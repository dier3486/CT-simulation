function [dataflow, prmflow, status] = reconnode_error(dataflow, prmflow, status)
% support node, to active an error
% [dataflow, prmflow, status] = reconnode_error(dataflow, prmflow, status);

% parameters set in pipe
errorprm = prmflow.pipe.(status.nodename);

status.jobdone = false;
if isfield(errorprm, 'errorcode')
    status.errorcode = errorprm.errorcode;
else
    status.errorcode = 1;
end
if isfield(errorprm, 'errormsg')
    status.errormsg = errorprm.errormsg;
else
    status.errormsg = [];
end

end
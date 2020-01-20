function [dataflow, prmflow, status] = reconnode_datasave(dataflow, prmflow, status)
% support node, save data (to .mat)
% [dataflow, prmflow, status] = reconnode_datasave(dataflow, prmflow, status);

% parameters set in pipe
backupprm = prmflow.pipe.(status.nodename);

% TBC


% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
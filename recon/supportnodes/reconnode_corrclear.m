function [dataflow, prmflow, status] = reconnode_corrclear(dataflow, prmflow, status)
% support node, clear corr tables
% [dataflow, prmflow, status] = reconnode_corrclear(dataflow, prmflow, status);

prmflow.corrtable = struct();

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
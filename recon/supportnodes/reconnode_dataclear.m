function [dataflow, prmflow, status] = reconnode_dataclear(dataflow, prmflow, status)
% support node, clear
% [dataflow, prmflow, status] = reconnode_dataclear(dataflow, prmflow, status);

dataflow = struct();
prmflow = struct();
status = struct();

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
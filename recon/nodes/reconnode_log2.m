function [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status)
% recon node, log2
% [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status);

% Z0
if isfield(prmflow.system, 'DBBzero')
    Z0 = prmflow.system.DBBzero;
else
    Z0 = 16384;
end
% log2
dataflow.rawdata = dataflow.rawdata - Z0;
dataflow.rawdata(dataflow.rawdata<=0) = nan;
dataflow.rawdata = -log2(dataflow.rawdata) + log2(single(dataflow.rawhead.Integration_Time));

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
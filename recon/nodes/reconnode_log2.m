function [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status)
% recon node, log2
% [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status);

% Z0
if isfield(prmflow, 'system') && isfield(prmflow.system, 'DBBzero')
    Z0 = prmflow.system.DBBzero;
else
    Z0 = 16384;
end
% offset
if isfield(dataflow, 'offset')
    dataflow.rawdata = dataflow.rawdata - mean(dataflow.offset.rawdata, 2);
else
    dataflow.rawdata = dataflow.rawdata - Z0;
end
dataflow.rawdata(dataflow.rawdata<=0) = nan;
% log2
dataflow.rawdata = -log2(dataflow.rawdata) + log2(single(dataflow.rawhead.Integration_Time));

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
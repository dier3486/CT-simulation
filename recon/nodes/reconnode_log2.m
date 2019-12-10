function [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status)
% recon node, log2
% [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status);

% log2
Z0 = 16384;
rawdata = rawdata - Z0;
rawdata(rawdata<=0) = nan;
rawdata = -log2(rawdata) + log2(single(rawhead.Integration_Time));


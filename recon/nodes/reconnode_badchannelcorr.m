function [dataflow, prmflow, status] = reconnode_badchannelcorr(dataflow, prmflow, status)
% recon node, log2
% [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status);

% input bad channel index

% calibration table
if isfield(prmflow.corrtable, status.nodename)
    badcorr = prmflow.corrtable.(status.nodename);
end




% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
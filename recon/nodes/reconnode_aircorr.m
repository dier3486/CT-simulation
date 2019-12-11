function [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status)
% recon node, air correction
% [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status);

% calibration table
aircorr = prmflow.corrtable.(status.nodename);



1;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
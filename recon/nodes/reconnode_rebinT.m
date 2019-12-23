function [dataflow, prmflow, status] = reconnode_rebinT(dataflow, prmflow, status)
% recon node, transpose of the rebin
% [dataflow, prmflow, status] = reconnode_rebinT(dataflow, prmflow, status);
% use to find out the inverse of the rebin

% rebin prepare
[prmflow, ~] = reconnode_rebinprepare(prmflow, status);

% TBC
% % Azi rebin
% [dataflow, prmflow, ~] = reconnode_Azirebin(dataflow, prmflow, status);
% 
% % Radial rebin
% [dataflow, prmflow, ~] = reconnode_Radialrebin(dataflow, prmflow, status);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
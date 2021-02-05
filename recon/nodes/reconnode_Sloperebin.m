function [dataflow, prmflow, status] = reconnode_Sloperebin(dataflow, prmflow, status)
% recon node, Axial 'slope' rebin 
% [dataflow, prmflow, status] = reconnode_Sloperebin(dataflow, prmflow, status);
% Support (X)DFS, gantry tilt, NO QDO

%1 slope rebin prepare
[prmflow, ~] = reconnode_sloperebinprepare(prmflow, status);

%2 slope fan-Radial
[dataflow, prmflow, ~] = reconnode_SlopeRadialrebin(dataflow, prmflow, status);

%3 slope Azi
[dataflow, prmflow, ~] = reconnode_SlopeAzirebin(dataflow, prmflow, status);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
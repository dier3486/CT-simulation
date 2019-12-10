% recon test script for Axial

% clear;
% path
addpath(genpath('../'));

reconxml = 'D:\matlab\data\simulation\recon_series1.xml';
reconcfg = readcfgfile(reconxml);
if ~iscell(reconcfg.recon)
    reconcfg.recon = {reconcfg.recon};
end

dataflow = struct();
prmflow = strcut();
status = struct();
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'initial');

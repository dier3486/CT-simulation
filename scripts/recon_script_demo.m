% recon demo script for Axial

% clear;
% path
% addpath(genpath('../'));

% load recon xml
% reconxml = 'D:\matlab\data\simulation\recon_series1.xml';
reconxml = 'E:\data\simulation\recon_series1.xml';
root = readcfgfile(reconxml);
if ~iscell(root.recon)
    root.recon = {root.recon};
end
% try 1st series
status.reconcfg = root.recon{1};
status.series_index = 1;

% initial
dataflow = struct();
prmflow = struct();
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'initial');

% load corr
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'loadcorrs');

% read raw data
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'loadrawdata');

% log2
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'log2');

% air
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Air');

% beamharden
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Beamharden');

% Housefield
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Housefield');

% rebin
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Axialrebin');

% FBP
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'FBP');


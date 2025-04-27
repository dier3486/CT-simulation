% recon demo script for Axial

% clear;
% path
% addpath(genpath('../'));

% load recon xml
% reconxml = 'D:\matlab\data\simulation\recon_series1.xml';
reconxml = 'E:\data\rawdata\bhtest\recon_series1.xml';
% reconxml = 'D:\data\simulation\aircali_test\recon_series1.xml';

root = readcfgfile(reconxml);
if ~iscell(root.recon)
    root.recon = num2cell(root.recon);
end
% try 1st series
status.reconcfg = root.recon{1};
status.seriesindex = 1;

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

% bad channel
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Badchannel');

% beamharden
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Beamharden');

% Hounsefield
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Hounsefield');

% rebin
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Axialrebin');

% prmflow.pipe.Filter.fillup = true;
% filter
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Filter');

% BP
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Backproject');

% % FBP
% [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'FBP');


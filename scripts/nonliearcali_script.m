% sample script for nonliear calibration

% clear;
% path
% addpath(genpath('../'));

% load recon xml
% reconxml = 'D:\matlab\data\simulation\recon_series1.xml';
reconxml = 'E:\data\rawdata\bhtest\cali_nonlinear1.xml';
% reconxml = 'E:\data\simulation\recon_series1.xml';

root = readcfgfile(reconxml);
if ~iscell(root.recon)
    root.recon = {root.recon};
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

% Housefield
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Housefield');

% record rawdata
% tmp
dataflow.rawdata_bh = dataflow.rawdata;
dataflow.rawdata_bh = reshape(dataflow.rawdata, prmflow.recon.Npixel, prmflow.recon.Nslice, prmflow.recon.Nview);

% rebin
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Axialrebin');

% ideal projection
prmflow.pipe.watergoback = struct();
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'watergoback');

% inverse rebin
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'Inverserebin');





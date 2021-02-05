% 'slope' rebin test code
% Atilt transform

load('E:\data\simulation\TM\test\tiltrb_test1.mat');
prmflow.recon.gantrytilt = prmflow.protocol.gantrytilt*(pi/180);

gpuDevice;

[dataflow, prmflow, status] = reconnode_Sloperebin(dataflow, prmflow, status);

% Filter
status.nodename = 'Filter';
[dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status);

% BP prepare
status.nodename = 'Backprojection';
[dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status);



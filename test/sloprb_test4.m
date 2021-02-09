% 'slope' rebin test code
% Atilt transform

load('E:\data\simulation\TM\test\tiltrb_test1.mat');
% fillup missed values
prmflow.recon.gantrytilt = prmflow.protocol.gantrytilt*(pi/180);
% prmflow.protocol.reconFOV = 300;
prmflow.protocol.imagesize = 768;
prmflow.pipe.Backprojection.method = '3D';
prmflow.pipe.Backprojection.FOV = 300;

gpuDevice;
% rebin
[dataflow, prmflow, status] = reconnode_Sloperebin(dataflow, prmflow, status);

% Filter
status.nodename = 'Filter';
[dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status);

% BP prepare
status.nodename = 'Backprojection';
[dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status);

% BP
[dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status);



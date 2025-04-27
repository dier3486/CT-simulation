function [dataflow, prmflow, status] = reconnode_multiaxialmeanprepare(dataflow, prmflow, status)
% prepare of mean multi-rotation of rawdata
% [dataflow, prmflow, status] = reconnode_multiaxialmeanprepare(dataflow, prmflow, status);

% parameters set in pipe
nodeprm = prmflow.pipe.(status.nodename);

% if to mean all the shots
if isfield(nodeprm, 'meanshots')
    meanshots = nodeprm.meanshots;
else
    meanshots = false;
end

% parameters to use in prmflow
if meanshots
    Nshot = 1;
else
    Nshot = prmflow.raw.Nshot;
end
Nview = prmflow.raw.Nview;
Nviewprot = prmflow.raw.Nviewprot;
prmflow.raw.Nmulti = Nview/Nviewprot;

% prmflow.raw after mean ??
prmflow.raw.Nview = Nviewprot * Nshot;
prmflow.raw.Nshot = Nshot;
prmflow.raw.viewpershot = Nviewprot;
prmflow.raw.viewnumber = prmflow.raw.viewnumber/prmflow.raw.Nmulti;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
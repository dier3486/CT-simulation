function [dataflow, prmflow, status] = reconnode_rotationmeanprepare(dataflow, prmflow, status)
% prepare node, prepare of rotation mean
% [dataflow, prmflow, status] = reconnode_rotationmeanprepare(dataflow, prmflow, status);

% reset parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

prmflow.pipe.(nodename).nominalnumber = prmflow.raw.Nviewprot;
prmflow.pipe.(nodename).viewskip = 0;
if ~isfield(nodeprm, 'overshots') || isempty(nodeprm.overshots)
    prmflow.pipe.(nodename).overshots = true;
end

% call datameanprepare
[dataflow, prmflow, status] = reconnode_datameanprepare(dataflow, prmflow, status);


end
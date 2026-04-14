function [dataflow, prmflow, status] = reconnode_pipelinestuckprepare(dataflow, prmflow, status)
% prepare node, pipeline stuck prepare
% [dataflow, prmflow, status] = reconnode_pipelinestuckprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% pipeline-stuck is a hard coded reconnode_databuffer
% to hard code these values
if ~isfield(nodeprm, 'alldata') && ~isfield(nodeprm, 'bufferfields')
    prmflow.pipe.(nodename).alldata = true;
end

prmflow.pipe.(nodename).copytodataflow = true;
prmflow.pipe.(nodename).stuck = true;

% and trun to databuffer-prepare
[dataflow, prmflow, status] = reconnode_databufferprepare(dataflow, prmflow, status);

end
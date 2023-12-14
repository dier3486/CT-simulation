function [dataflow, prmflow, status] = reconnode_pipelinestuck(dataflow, prmflow, status)
% support node, to link a pipeline node with a non-pipeline node
% [dataflow, prmflow, status] = reconnode_pipelinestuck(dataflow, prmflow, status);

% pipeline stuck is a hard coded reconnode_databuffer, see reconnode_pipelinestuckprepare

% call reconnode_databuffer
[dataflow, prmflow, status] = reconnode_databuffer(dataflow, prmflow, status);

end
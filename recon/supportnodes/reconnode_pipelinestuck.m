function [dataflow, prmflow, status] = reconnode_pipelinestuck(dataflow, prmflow, status)
% support node, to link a pipeline node with a non-pipeline node
% [dataflow, prmflow, status] = reconnode_pipelinestuck(dataflow, prmflow, status);
% pipeline stuck is a hard coded reconnode_databuffer, see reconnode_pipelinestuckprepare

% not prepared?
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_pipelinestuckprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% call reconnode_databuffer
[dataflow, prmflow, status] = reconnode_databuffer(dataflow, prmflow, status);

end
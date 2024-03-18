function [dataflow, prmflow, status] = reconnode_null(dataflow, prmflow, status)
% support node, not exist null node
% [dataflow, prmflow, status] = reconnode_null(dataflow, prmflow, status);

error('A NULL node does not exist! It can not be configured or called by any way!');

end
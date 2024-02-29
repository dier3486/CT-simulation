function [privpool, privdata] = readdatafrominput(privpool, privdata, inputpool, inputdata)
% pipe-line function to read data from an inputpool to a private pool

% e.g.
% privpool = dataflow.buffer.(nodename).outpool;
% privdata = dataflow.buffer.(nodename).outpool.data;
% inputpool = status.pipepool.(nodename);
% inputdata = dataflow.pipepool.(nodename);

% datasize
if isfield(inputpool, 'AvailPoint')
    ReadendPoint = inputpool.AvailPoint;
else
    ReadendPoint = inputpool.WritePoint;
end
datasize_toread = max(0, ReadendPoint - inputpool.ReadPoint);

% check the buffer left
privpoolsize = privpool.poolsize - privpool.WritePoint + 1;
datasize_toread = min(datasize_toread, privpoolsize);

if isfield(privpool, 'datafields')
    datafields = privpool.datafields;
else
    datafields = {};
end
% copy pool to buffer.outpool
[privdata, ~] = pooldatacopy(inputdata, privdata, inputpool.ReadPoint, privpool.WritePoint, datasize_toread, {}, true);

    % move current pool's read point
    status.pipepool.(nodename).ReadPoint = status.pipepool.(nodename).ReadPoint + datasize_toread; 
    % move the buffer.outpool's write point
    dataflow.buffer.(nodename).WritePoint = dataflow.buffer.(nodename).WritePoint + datasize_toread;
    % Note: the buffer.(nodename).outpool is a private buffer, not the next pool. The data in buffer.(nodename).outpool will be
    % copied to the next pool after the correcion works done.


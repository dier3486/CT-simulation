function [privpool, privdata, inputpool, readnumber] = readdatafrominput(privpool, privdata, inputpool, inputdata, datasize_toread)
% pipe-line function to read data from an inputpool to a private pool
% [privpool, privdata, inputpool, readnumber] = readdatafrominput(privpool, privdata, inputpool, inputdata);
% e.g.
%   privpool = dataflow.buffer.(nodename).outpool;
%   privdata = dataflow.buffer.(nodename).outpool.data;
%   inputpool = status.pipepool.(nodename);
%   inputdata = dataflow.pipepool.(nodename);
% the readnumber is the view number was read.

% skip if the privpool is empty
if isempty(privpool)
    readnumber = 0;
    return;
end

% datasize to read
if isfield(inputpool, 'AvailPoint')
    ReadendPoint = inputpool.AvailPoint;
else
    ReadendPoint = inputpool.WritePoint;
end
if nargin < 5
    datasize_toread = max(0, ReadendPoint - inputpool.ReadPoint);
else
    datasize_toread = max(min(datasize_toread, ReadendPoint - inputpool.ReadPoint), 0);
end

% WriteEnd
if isfield(privpool, 'WriteEnd')
    WriteEnd = min(privpool.WriteEnd, privpool.poolsize);
else
    WriteEnd = privpool.poolsize;
end

% check the buffer left
privpoolsize = WriteEnd - privpool.WritePoint + 1;
datasize_toread = min(datasize_toread, privpoolsize);

% accept fields
if isfield(privpool, 'datafields')
    datafields = privpool.datafields;
else
    % default is all
    datafields = {};
end

% copy pool to buffer.outpool
[privdata, readnumber] = ...
    pooldatacopy(inputdata, privdata, inputpool.ReadPoint, privpool.WritePoint, datasize_toread, datafields, true);

% move input pool's read point
inputpool.ReadPoint = inputpool.ReadPoint + readnumber;
if isfield(inputpool, 'ReadViewindex')
    inputpool.ReadViewindex = inputpool.ReadViewindex + readnumber;
end

% move the buffer.outpool's write point
privpool.WritePoint = privpool.WritePoint + readnumber;

end


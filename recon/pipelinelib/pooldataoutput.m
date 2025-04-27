function [currpool, nextpool, Rinfo] = pooldataoutput(currpool, nextpool, viewcommon, maxcopynum)
% to copy the data from currpool to nextpool
% a simplified version in using pooldatacopy
% [currpool, nextpool, Rinfo] = pooldataoutput(currpool, nextpool, viewcommon, maxcopynum);
% or, [currpool, nextpool, Rinfo] = pooldataoutput(currpool, nextpool);

% the view number must be integer times to viewcommon
if nargin < 3 || isempty(viewcommon)
    viewcommon = 1;
end

if nargin < 4
    maxcopynum = inf;
end

% datasize to write is the data size can be read in currpool
datasize_towrite = poolreadable(currpool, viewcommon);
datasize_towrite = min(datasize_towrite, maxcopynum);
% which will be returned in Rinfo
Rinfo.datasize_towrite = datasize_towrite;

% datasize can write is the buffer size left in nextpool
datasize_canwrite = poolpspaceleft(nextpool, viewcommon);

if ~isempty(nextpool)
    Rinfo.writenum = min(datasize_towrite, datasize_canwrite);
    Rinfo.origWritePoint = nextpool.WritePoint;
else
    % when the nextpool is NULL, the copy will play like usual
    Rinfo.writenum = datasize_towrite;
    Rinfo.origWritePoint = 0;
    % move the points
    [currpool, nextpool] = movepointsaftercopy(currpool, nextpool, Rinfo.writenum);
    return;
end

% pool data copy
if isfield(nextpool, 'dynamicfields')
    dynamicfields = nextpool.dynamicfields;
else
    dynamicfields = false;
end
nextpool.data = pooldatacopy(currpool, currpool.data, nextpool, nextpool.data, Rinfo.writenum, currpool.datafields, ...
    dynamicfields);
% move the points
[currpool, nextpool] = movepointsaftercopy(currpool, nextpool, Rinfo.writenum);

% WriteEnd to WriteEnd of nextpool
WriteClosing = nextpool.WriteEnd - nextpool.WritePoint;

% read/write lock for circulatemode
if nextpool.circulatemode && isfield(nextpool, 'WriteStuck')
    nextpool.WriteStuck = (Rinfo.writenum>0) & (WriteClosing==1);
end
% we don't lock the reading

% flag of jobdone
if (Rinfo.datasize_towrite > 0) && (Rinfo.datasize_towrite > Rinfo.writenum)
    % nextpool could be stucked
    Rinfo.jobflag = 2;
elseif Rinfo.datasize_towrite > 0 && (Rinfo.datasize_towrite == Rinfo.writenum)
    % wrote something
    Rinfo.jobflag = 1;
elseif Rinfo.datasize_towrite == 0
    % nothing wrote
    Rinfo.jobflag = 3;
else
    % totally stucked
    Rinfo.jobflag = 6;
end

end


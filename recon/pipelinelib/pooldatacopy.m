function [nextdata, writenum] = pooldatacopy(currpool, currdata, nextpool, nextdata, writenum, ...
    copyfields, force_flag, method)
% to copy the data from currunt data to next data in pipeline, circulated patch
%   nextdata = pooldatacopy(currpool, currdata, nextpool, nextdata, writenum, copyfields, force_flag);
% or
%   nextdata = pooldatacopy(currpool, currdata, nextpool, nextdata, writenum);
% It will be moved to pooldatacopy as an inner subfucntion

if nargin < 7
    force_flag = false;
end

if nargin < 6 || isempty(copyfields)
    if force_flag
        copyfields = currpool.datafields;
    else
        copyfields = nextpool.datafields;
    end
end

if nargin < 8
    method = 'overwrite';
end

% the reading data's index
readindex = currpool.ReadPoint : currpool.ReadPoint + writenum - 1;
if currpool.circulatemode
    readindex = mod(readindex - 1, currpool.poolsize) + 1;
    % the writenum shall <= currpool.poolsize
end
% the writing data' index
writeindex = nextpool.WritePoint : nextpool.WritePoint + writenum - 1;
if nextpool.circulatemode
    writeindex = mod(writeindex - 1, nextpool.poolsize) + 1;
    % the writenum shall <= nextpool.poolsize
end

% the out range index will be ignore
s = (writeindex < 1) | (writeindex > nextpool.poolsize);
if any(s)
    writeindex(s) = [];
    readindex(s) = [];
end
s = (readindex < 1) | (readindex > currpool.poolsize);
if any(s)
    writeindex(s) = [];
    readindex(s) = [];
end
writenum = length(writeindex);
1;
nextdata = poolhardcopy(nextdata, currdata, writeindex, readindex, copyfields, method);

end